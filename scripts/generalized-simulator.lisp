;;;; generalized-simulator.lisp --- Implement the MSCR-JC(k) process on
;;;;                                arbitrary binary trees. Currently under
;;;;                                construction.
;;
;; Author: max hill 
;; (Last updated 2022-05-24)

;; DESCRIPTION: Here we attempt to implement a simulator which takes as input a
;; binary tree T (ideally in newick tree format, with branch lengths) and output
;; an MSA in the form of an alignment file generated according to the MSCR
;; process on T.
;;
;; At present, both the JC69 process and the MSCR process have been implemented
;; for arbitrary trees. These are two ends of the pipeline, but they are not
;; connected yet. Inference is also not implemented.
;;
;; Remains to do:
;;
;; 1. (low priority) Implement a newick -> lisp tree converter
;;
;; 2. (high priority) Write a program to generate marginal gene trees from the
;;    output of mscr which can then be used as input for the mutation simulator.
;;    One option: to generate a marginal gene tree for each of the k sites, we
;;    can run (compute-marginal-tmrcas i j output-edges k) for all distinct
;;    pairs of leaves i and j. Then we'll have a big matrix of coalescent times
;;    and we can just lookup the times to create the marginal gene trees. Or
;;    better, we can write a new function that takes as input the set
;;    output-edges, a site i∈{1,...k}, and outputs a gene tree in the
;;    appropriate form for evolve-down-tree. Then run evolve down tree and save
;;    the result as a column in the MSA.
;;
;; 3. (high priority) Implement the inference methods on the new MSAs. In
;;    particular, we want to implement the four-point method. We could also try
;;    estimating branch lengths, as this would be a key test of one of our
;;    paper's conclusions.
;;
;; 4. (low priority) address performance issues arising from use of integer
;;    intervals.
;;
;; Sources: The following resource was extremely helpful:
;; https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html
;;
;;______________________________________________________________________________
;;
;; Part 1. Constructors
;;______________________________________________________________________________
;;
;; There are two kinds of binary trees: (1) nodes and (2) leaves. 
;;
;; 1. A *leaf* is a list with two components: a string (i.e. the label) followed
;;    by a list of 4 elements (1: a number representing the distance from the
;;    parent, 2: a number representing the time of the leaf node, 3: a number,
;;    representing the mutation rate on that edge, and 3: a number, representing
;;    the recombination rate on that edge). For example, ("a" (2.3 .2 0 .1 1.2))
;;
;;    ("label" (:edge-length :vertex-age :mutation-rate :recomb-rate))
;;
;; 2. A *node* is a list with three components: a leaf, a left subtree and a
;;    right subtree: (:leaf :left-subtree :right-subtree).
;;
;; The following discussion assumes that the leaves of the tree are at the
;; bottom and the root is at the top.
;;
;; - :recomb-rate is the rate of recombination along the edge extending upward
;;   from the tree's top vertex.
;;
;; - :mutation-rate is the mutation rate along the edge extending upward from
;;   the tree's top vertex
;;
;; - :dis-from-parent is the length of the edge connecting the vertex to its
;;   parent.
;;
;; - :vertex-age is the age of the vertex


;; A note on the terminology used in this simulator. First, the phylogenetic
;; trees (i.e. the species tree) here are binary and are always thought to be
;; oriented so that the root is at the top and the leaves are at the bottom. The
;; present time is t=0, and so more recent things are located lower on the tree,
;; and "time" really means "age". Second, the phylogenetic tree is thought of as
;; a structured collection of populations; in particular, each vertex and the
;; edge directly above it is regarded as a "population". The "start time" of a
;; population is the age of the vertex. The end time of a population is the
;; start time plus the length of the edge. Here we regard the edge extending up
;; from the vertex to be a 'population.' Hence, the population-start-time is the
;; age of the vertex, and the population-end-time is the age of its
;; parent.Third, mutation and recombination rates are associated individually
;; with each population -- they need not all be the same everywhere on the tree.
;; Fourth, there are two kinds of populations: "leaf" populations, which consist
;; only of a list containing an plist of parameters, and "node" populations,
;; which correspond to internal vertices of the phylogenetic tree. Node
;; populations consist of a list containing three things: an plist of
;; parameters, a left subtree, and a right subtree (which are both themselves
;; either leaf or node populations.

(defun make-leaf (leaf-label distance-from-parent mutation-rate recombination-rate)
  "Create a leaf."
  (list
   (list :population-label leaf-label
         :dist-from-parent distance-from-parent
         :mutation-rate mutation-rate
         :recomb-rate recombination-rate
         :parent-label nil
         :population-start-time nil
         :population-end-time nil
         :numeric-label nil)))

(defun make-node (node-label distance-from-parent mutation-rate
                  recombination-rate left-subtree right-subtree)
  "Create a node."
  (list
   (list :population-label node-label
         :dist-from-parent distance-from-parent
         :mutation-rate mutation-rate
         :recomb-rate recombination-rate
         :parent-label nil
         :population-start-time nil
         :population-end-time nil
         :numeric-label nil)
   left-subtree
   right-subtree))

;;______________________________________________________________________________
;;
;; Part 2. Recognizers
;;______________________________________________________________________________

(defun leafp (tree)
  "Return t if tree is a leaf and nil otherwise."
  (and (listp tree)
       (= 1 (length tree))
       (listp (first tree))))

(defun nodep (tree)
  "Return t if tree is a node and nil otherwise."
  (and (listp tree)
       (= 3 (length tree))
       (listp (first tree))
       (or (nodep (second tree))
           (leafp (second tree)))
       (or (nodep (third tree))
           (leafp (third tree)))))

;; The next function tests whether a population is the root
(defun rootp (tree)
  "Test whether the top population of tree is the root population."
  (not (get-parameter :parent-label tree)))


;;______________________________________________________________________________
;;
;; Part 2. Selectors
;;______________________________________________________________________________

(defun get-parameter (parameter population)
  "Retrieve specific parameters from a population"
  (getf (first population) parameter))

(defun get-parameter-plist (population)
  "Retrieve the property list of parameters for the given population."
  (first population))

(defun set-population-parameter (parameter value population)
  "Set an already-existing population parameter to given value."
  (setf (getf (first population) parameter) value))

(defun add-population-parameter (parameter-key parameter-value tree)
  "Adds parameter with given value to the plist of a single population. Does
not affect subtrees."
  (let* ((updated-parameter-plist
           (list* parameter-key
                  parameter-value
                  (get-parameter-plist tree))))
    (if (leafp tree)
        (list updated-parameter-plist)
        (list updated-parameter-plist
              (left-subtree tree)
              (right-subtree tree)))))

(defun left-subtree (tree)
  "Return the left subtree of tree. If tree is a leaf, return nil."
  (if (leafp tree) nil (second tree)))

(defun right-subtree (tree)
  "Return the right subtree of tree. If tree is a leaf, return nil."
  (if (leafp tree) nil (third tree)))

(defun get-population-name (tree)
  "Return the name (i.e. label) of the population. Works for both leaves and internal
nodes."
  (get-parameter :population-label tree))

(defun get-parent (tree child-name)
  "Returns the parent node of the vertex with label child-name, which is a
string. Usage: (get-parent *n3* ``ab''). To understand this recursive function,
one needs to understand two properties of cond: first, cond terminates when the
first condition (evaluated in order) evaluates to a non-nil result. Second, if
test1 is nonnil, then a cond condition of the form ((test1)) returns the value
of test1."
  (cond ((leafp tree) "leafy")
        ((or (equal (get-population-name (left-subtree tree)) child-name)
             (equal (get-population-name (right-subtree tree)) child-name))
         tree)
        ((get-parent (left-subtree tree) child-name))
        ((get-parent (right-subtree tree) child-name))))

(retrieve left-daugther population)
(retrieve right-daugther population)
(retrieve population label)
(retrieve population mutation-rate)
(retrieve population recomb-rate)
(retrieve population start-time)
(retrieve population end-time)
(retrieve population edge-length)
(retrieve parent population)


          




;;______________________________________________________________________________
;;
;; Part 3. Augmenting tree Parameters
;;______________________________________________________________________________

;; This section defines a function which compute additional information about
;; the tree and adds it to the parameter list.


;; First we define some counting functions to obtain the heighet as well as the
;; number of leaves on the binary tree.
(defun count-number-of-leaves (tree)
  "Count the number of leaves on a binary tree."
  (if (leafp tree)
      1
      (+ (count-number-of-leaves (left-subtree tree))
         (count-number-of-leaves (right-subtree tree)))))

(defun get-tree-height (tree)
  "Return the total height of the tree, i.e. from the top of the root edge to
  the maximally distant leaf."
  (get-tree-height-aux tree 0))

(defun get-tree-height-aux (tree x)
  "Auxilliary function for get-tree-height. The variable x is a dummy variable
  for adding up the total height of the tree."
  (let ((dist-from-parent (get-parameter :dist-from-parent tree)))
    (if (leafp tree)
        (+ x dist-from-parent)
        (max
         (get-tree-height-aux (right-subtree tree) (+ x dist-from-parent))
         (get-tree-height-aux (left-subtree tree) (+ x dist-from-parent))))))


;; The next function augments the tree parameter plists by recording information
;; about parent populations (for easy retreival) and calculating population ages
;; from edge lengths (which will aslo be convenient later).

(defun add-ages-and-parents
    (tree &optional (parent-age (get-tree-height tree)) (parent-label nil))
  "Augment the population parameter plists of tree with additional information
about the parents and ages of the populations. The optional variables are used
for recursion and should not be entered by the user."
  (let* ((dist-from-parent (get-parameter :dist-from-parent tree))
         (current-population-start-time (- parent-age dist-from-parent))
         (current-population-label (get-parameter :population-label tree)))
    (set-population-parameter :parent-label parent-label tree)
    (set-population-parameter :population-start-time current-population-start-time tree)
    (set-population-parameter :population-end-time parent-age tree)
    (unless (leafp tree)
      (add-ages-and-parents (left-subtree tree)
                            current-population-start-time
                            current-population-label)
      (add-ages-and-parents (right-subtree tree)
                            current-population-start-time
                            current-population-label))))

;; The next two functions augment the input tree parameter lists with a
;; structured population labeling system.

(defun add-numeric-labels (tree)
  "Add numeric labels to the tree parameter plists so that leaves are labeled
with the numbers 1,..,n, and the label of each other nodes is a list consisting
of the labels of its two children. This function both changes the tree and returns
the value of the updated tree. The changes made to the tree are idempotent."
  (progn
    (defparameter *leaf-number-tracker* 0)
    (add-numeric-labels-aux tree)
    (makunbound '*leaf-number-tracker*)
    tree))

(defun add-numeric-labels-aux (tree)
  "Recursive auxillary function for add-numeric-labels. Note: incf not only
increments the variable, but also outputs the new value."
  (set-population-parameter
   :numeric-label
   (if (leafp tree)
       (incf *leaf-number-tracker*)
       (list (add-numeric-labels-aux (left-subtree tree))
             (add-numeric-labels-aux (right-subtree tree))))
   tree))


;; Finally, the key function of this section is the following.

(defun augment-tree-parameters (tree)
  "The input should be a tree constructed with the functions 'make-node' and
'make-leaf'. This function is idempotent."
  (progn (add-ages-and-parents tree)
         (add-numeric-labels tree)))


;;______________________________________________________________________________
;;
;; Part 4. Tree examples for testing
;;______________________________________________________________________________
;; From here on we will assume that all trees have vertex times, i.e. were
;; processed with the function 'add-age-parameters-to-tree'

;; Example 1. A leaf
(defparameter *l1
  (make-leaf "A" .1 1 4))
(defparameter *lv1 (augment-tree-parameters *l1))
(add-numeric-labels *lv1)


;; Example 2. A node (unbalanced quartet):
(defparameter *n1
  (make-node "root" 0 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv1 (augment-tree-parameters *n1))
(add-numeric-labels *nv1)


;; Example 3. A node (unbalanced quartet):
(defparameter *n2
  (make-node "root" 999 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv2 (augment-tree-parameters *n2))
(add-numeric-labels *nv2)


;; Example 4. A node (balanced quartet):
(defparameter *n3
  (make-node "root" 999 .1 .2
             (make-node "ab" 2 .1 .2
                        (make-leaf "a" 2 .1 .2)
                        (make-leaf "b" 2 .1 .2))
             (make-node "cd" 1 .1 .2
                        (make-leaf "c" 3 .1 .2)
                        (make-leaf "d" 2 .1 .2))))
(defparameter *nv3 (augment-tree-parameters *n3))
(add-numeric-labels *nv3)


;; Example 5. A 6-taxa tree
(defparameter *bigbad
    (make-node "root" 1 2 3
               (make-node "ABC" 1 2 3 
                          (make-node "AB" 1 2 3
                                     (make-leaf "A" 1 2 3)
                                     (make-leaf "B" 1 2 3))
                          (make-leaf "C" 1 2 3))
               (make-node "D123" 1 2 3
                          (make-node "D12" 1 2 3
                                     (make-leaf "D1" 1 2 3)
                                     (make-leaf "D2" 1 2 3))
                          (make-leaf "D3" 1 2 3))))

(defparameter *bigbad (augment-tree-parameters *bigbad))
(add-numeric-labels *bigbad)


;;______________________________________________________________________________
;;
;; Part 6. Mutation Simulator.
;;______________________________________________________________________________
;; This part simulates the JC process on a (general) tree.
;;
;; For simplicity, we will implement nucleotides here as the numbers 0,1,2,3
;; rather than A,T,C,G. Everything will assume JC69 process.
;;
;; Each edge of the tree has an associated mutation rate. This is the rate at
;; which a nucleotide leaves its current state according to the JC69 process.

;; The next function is used to generate a nucleotide at the root of the tree.

(defun draw-random-nucleotide ()
  "Return a random number from the set {0,1,2,3}."
  (random 4))

;; The next (highly inefficient) function creates an interval of integers. Its
;; performance is very poor when k is large. With a little bit of thought and
;; effort, I could addres this issue. (low priority)

(defun interval (a b)
  "Construct an 'interval' of integers from a to b, inclusive. Note that a and b
must both be integers."
  (loop for i from a to b collecting i))

;; We'll need a way to perform nucleotide substitutions when they are determined
;; to occur. This functionality is provided by the next function.

(defun convert-nucleotide-to-letter (number)
  (cond ((= number 0) "A")
        ((= number 1) "T")
        ((= number 2) "C")
        ((= number 3) "G")))
        

(defun implement-substitution (current-nucleotide &optional population-name)
  "Implements a substitution to a *different* nucleotide. Specifically, return a
random number from the set {0,1,2,3}\{current-nucleotide}. It does this by
adding 1d3 to the current nucleotide and then reducing the result modulo 4."
  (declare (integer current-nucleotide))
  (let ((new-nucleotide (mod (+ current-nucleotide 1 (random 3)) 4)))
    (format t "~%Substitution from ~a to ~a in population ~a~%"
            current-nucleotide
            new-nucleotide
            population-name)
    new-nucleotide))

;; Substitutions will be performed on each edge with some probability, which is
;; what the next function does. It also contains some calls to format; these are
;; currently included for testing purposes only.

(defun evolve-down-edge (tree parent-nucleotide)
  "Implements substitutions along an edge, according to JC69 process. If no
substitutions, then return the parent nucleotide. Assumes JC procees as
specified in Dasarathy, Mossel, Nowak, Roch 'Coalescent-based species tree
estimation: a stochastic Farris transform'. Output takes the form of a
pair (label-name nucleotide)."
  (let* ((population-name (get-population-name tree))
         (edge-length (get-parameter :dist-from-parent tree))
         (mutation-rate (get-parameter :mutation-rate tree))
         (substitution-probability
           (* .75 (- 1 (exp (* (/ -4 3) mutation-rate edge-length)))))
         (new-nucleotide
           (if (< (random 1d0) substitution-probability)
               (implement-substitution parent-nucleotide population-name)
               (progn (format t "~%No substitution in population ~a~%"
                              population-name) parent-nucleotide))))
    (progn
      (if (leafp tree)
          (format t "leaf ~a: ~a ~%" population-name new-nucleotide)
          (unless (rootp tree)
            (format t "ancestor ~a: ~a ~%" population-name new-nucleotide)))
      (list population-name new-nucleotide))))

;; The next two functions are the workhorses of our sequence simulator. They use
;; the above functions to recursively implement the JC69 process on the input
;; tree.


(defun evolve-down-tree-aux (tree parent-nucleotide x)
  "Recursive auxilliary function for imlement-jc-process. The variable x is the
accumulator variable use to construct the output list."
  (if (leafp tree)
      (cons (evolve-down-edge tree parent-nucleotide) x)
      (let ((new-parent-nucleotide
             (second (evolve-down-edge tree parent-nucleotide))))
        (evolve-down-tree-aux (left-subtree tree)
                              new-parent-nucleotide
                              (evolve-down-tree-aux (right-subtree tree)
                                                    new-parent-nucleotide
                                                    x)))))

(defun implement-jukes-cantor-process (tree)
  "Recursive implmentation of JC69 on the tree. Starts with random nucleotide at
the tip of the edge above the root. To make this start exactly at the root (by
which I mean the mrca of the samples) you should set the 'dist-from-parent' in
the root vertex to zero."
  (let ((root-state (draw-random-nucleotide)))
    (progn (format t "Root state: ~a ~%" root-state)
           (evolve-down-tree-aux tree root-state nil))))

;; It works! Example usage:
;; (implement-jukes-cantor-process *n1)

;; The following function allows us to avoid the strange backticks like those
;; used in our original simulator, simulate-three-species

(defun make-leaf-sample (i n-total-leaf-number k-sequence-length leaf)
  "makes the i-th initial sample for a phylogenetic tree with a number of taxa
equal to n-total-leaf-number, and when sequences are length k-sequence-length"
  (let ((new-lineage
          (list
           (list
            (cons (get-parameter :population-start-time leaf)
                  (append (make-list (- i 1))
                          (cons (interval 1 k-sequence-length)
                                (make-list (- n-total-leaf-number i)))))))))
    (progn (format t "~%~%Sample ~a created in population ~a at time ~a~%One lineage added: ~a"
                   i
                   (get-population-name leaf)
                   (get-parameter :population-start-time leaf)
                   (first (first new-lineage)))
           new-lineage)))

(defun merge-output-edge-sets (child-1 child-2)
  "Combine the output edge sets child-1=(P_1,Q_1) and child-2=(P_2,Q_2) from two
daughter populations for entry into their parent population. Output a pair
of (P,Q) of active and inactive edges in the appropriate format for use in
arg-builder. Note that if one of the edge-sets is nil, this will just output the
other edge set. Recall that the set P consists of those lineages which are
active in the sense that they may still undergo coalescent or recombination
events. The seq Q contains lineages which have been deactivated as a result of
undergoing coalescent events -- but these can't just be thrown out because they
contain information about MRCAs."
  (list (union (first child-1) (first child-2))
        (union (second child-1) (second child-2))))

(defun get-number-of-active-lineages (edge-sets)
  "Determine how many lineages are 'active' in a population (a lineage is active
if it is a candidate for coalescence. Recall that when a lineage undergoes
recombination, it is rendered inactive, but two new active lineages are created.
In general, active lineages are stored in P and inactive lineages are stored in
Q. I can't remember off the top of my head but there might be an exception to
this rule.)"
  (length (first edge-sets)))

;; Load functions from old simulator
(load "/home/mutalisk/research/recombination-project/MSCR-simulator/scripts/simulator.lisp")


;; This is the auxillary function to the main simulator; it simulates the ARG
;; process in a single population. It is analogous to 'arg-builder' in the
;; previous work

(defun build-single-population-arg (ρ t_0 t_end population-name k-sequence-length
                                    edge-sets &optional (stop-at-mrca nil))
  "Build an ancestral recombination graph in a single population with fixed
start time t_0 and end time t_end."
  (let* ((number-of-active-lineages (get-number-of-active-lineages edge-sets))
	 (coales-rate (* .5 number-of-active-lineages
                         (- number-of-active-lineages 1)))
	 (recomb-rate (* .5 number-of-active-lineages ρ))
	 (total-rate (+ coales-rate recomb-rate))
	 (t_1 (+ t_0 (draw-exponential total-rate))))
    (if (or (> t_1 t_end) (and stop-at-mrca (= number-of-active-lineages 1)))
	edge-sets
	(build-single-population-arg
         ρ t_1 t_end
         population-name
         k-sequence-length
         (if (< (random 1d0) (/ recomb-rate total-rate))
	     (implement-recombination t_1 edge-sets k-sequence-length population-name)
	     (implement-coalescence t_1 edge-sets population-name))
	 stop-at-mrca))))


;; This is the main simulator. It is currently running without error, but the
;; output is not correct. It looks like some of the subroutines borrowed from
;; the earlier simulator need to be modified... Actually, I think I fixed those
;; issues now. I need to make mscr into an auxillary function mscr-aux, and then
;; write an intial function mscr which treats the root population differently
;; (i.e. by setting 'stop-at-mrca' to true when calling
;; 'build-single-population-arg'.
(defun start-mscr (species-tree k-sequence-length)
  (mscr (count-number-of-leaves species-tree)
        k-sequence-length
        species-tree
        nil))

(defun mscr (n k tree edge-sets)
  "Input comments: the variable 'n' is number of leaves on the initial (i.e.
full) tree, 'k' is the length of the sequences in base pairs. The input
'leaf-number-tracker' tracks the leaf number: its initial value is 1? and it
increments leaf-number-tracker by 1 each time a right subtree is recursively
evaluated. The initial value of edge-sets should be nil, I think."
(if (leafp tree)
    (make-leaf-sample (get-parameter :numeric-label tree) n k tree)
    (let* ((lineages-from-left-daughter (mscr n k (left-subtree tree) edge-sets))
           (lineages-from-right-daughter (mscr n k (right-subtree tree) edge-sets))
           (starting-lineages (merge-output-edge-sets lineages-from-left-daughter
                                                      lineages-from-right-daughter)))
      (build-single-population-arg
       (get-parameter :recomb-rate tree)
       (get-parameter :population-start-time tree)
       (get-parameter :population-end-time tree)
       (get-parameter :population-label tree) ; maybe use the numeric labels?
       k ; sequence length
       starting-lineages
       (rootp tree)))))


;; This is the main simulator. It is currently running without error, but the
;; output is not correct. It looks like some of the subroutines borrowed from
;; the earlier simulator need to be modified... Actually, I think I fixed those
;; issues now. I need to make mscr into an auxillary function mscr-aux, and then
;; write an intial function mscr which treats the root population differently
;; (i.e. by setting 'stop-at-mrca' to true when calling
;; 'build-single-population-arg'.
(defun mscr (species-tree k-sequence-length)
  (mscr-aux (count-number-of-leaves species-tree)
            k-sequence-length
            species-tree))

(defun mscr-aux (n k tree)
  "Recursive auxilliary function for mscr. The variable 'n' is number of leaves on the initial (i.e.
full) tree, 'k' is the length of the sequences in base pairs."
  (if (leafp tree)
      (make-leaf-sample (get-parameter :numeric-label tree) n k tree)
      (let* ((lineages-from-left-daughter (mscr-aux n k (left-subtree tree)))
             (lineages-from-right-daughter (mscr-aux n k (right-subtree tree)))
             (starting-lineages (merge-output-edge-sets lineages-from-left-daughter
                                                        lineages-from-right-daughter)))
        (build-single-population-arg (get-parameter :recomb-rate tree)
                                     (get-parameter :population-start-time tree)
                                     (get-parameter :population-end-time tree)
                                     (get-parameter :population-label tree) ; maybe use the numeric labels?
                                     k ; sequence length
                                     starting-lineages
                                     (rootp tree)))))


;; left off here :)


;; ((4.584407128516051d0 (1) (1) (1) (1))
;;  (4.497471958701756d0 NIL NIL (1) NIL)
;;  (4.025203326294506d0 (1) (1) NIL (1))
;;  (2.693733411688349d0 (1) (1) NIL NIL))

;; ((4.584407128516051d0 (2) (2) (2) (2))
;;  (4.497471958701756d0 NIL NIL (2) NIL)
;;  (4.025203326294506d0 (2) (2) NIL (2))
;;  (2.693733411688349d0 (2) (2) NIL NIL))

1: (0 2.7 4.0 4.6) 
2: (0 2.7 4.0 4.6)
4: (1 4.0 4.6)
3: (0 4.6)

;; The following function will be used to convert our time matrix into input for
;; the Jukes-Cantor process. 

(defun who-coalesced (time-matrix coal-time n-number-of-species)
  "Return the set of leaf labels whose ancestral lineages underwent coalesence
at the given coalescence time. Might be better named by something like 'common-ancestor' or 'population-name'"
  (loop for i from 1 to (1- n-number-of-species)
        appending (loop for j from (1+ i) to n-number-of-species
                        if (= (get-matrix-entry time-matrix i j) coal-time)
                          collect i into inner-output
                          and collect j into inner-output
                        finally (return inner-output))
          into outer-output
        finally (return (remove-duplicates outer-output))))

;; To use the above function, we will need a descending list of coalescent
;; times.

(defun get-list-of-coal-times (time-matrix n-number-of-species)
  "Return a sorted (descending) set of coalescent times for the n species whose
times are give by time-matrix."
  (loop for row from 1 to (1- n-number-of-species)
        appending (loop for column from (1+ row) to n-number-of-species
                        collecting (get-matrix-entry time-matrix row column)
                          into inner-output
                        finally (return inner-output))
          into outer-output
        finally (return (sort (remove-duplicates outer-output)
                              '>))))


(get-list-of-coal-times *time-matrix* 4)
(defparameter *ct* '(4.584407128516051d0 4.025203326294506d0 2.693733411688349d0))
(who-coalesced *time-matrix* (third (get-list-of-coal-times *time-matrix* 4)) 4)
(who-coalesced *time-matrix* (second (get-list-of-coal-times *time-matrix* 4)) 4)
(who-coalesced *time-matrix* (first (get-list-of-coal-times *time-matrix* 4)) 4)

(defun make-marginal-tree (descending-coal-times time-matrix n)
  (if (= 1 (length (who-coalesced *time-matrix* (first descending-coal-times) n)
      (make-leaf descending-coal-times 1 .1 0)
      (let* ((current-lineage
               (who-coalesced *time-matrix* (first descending-coal-times) n))
             (left-lineage
               (who-coalesced *time-matrix* (second descending-coal-times) n))
             (right-lineage
               (set-difference current-lineage right-lineage)))
        (make-node current-lineage 1 .1 0
                   (make-marginal-tree left-lineage time-matrix n)
                   (make-marginal-tree right-lineage time-matrix n)))))


;; Essentially, we want to make something like this:
;; (make-node 
;;  "1234" 0 .1 0 
;;  (make-node 
;;   "124" 
;;   (- 4.584407128516051d0 4.025203326294506d0)
;;   0.1 0 
;;   (make-node 
;;    "12" 
;;    (- 4.025203326294506d0 2.693733411688349d0)
;;    .1 0
;;    (make-leaf "1" 0 .1 0)
;;    (make-leaf "2" 0 .1 0))
;;   (make-leaf
;;    "4"
;;    (- 4.025203326294506d0 1)
;;    .1 0))
;;  (make-leaf
;;   "3"
;;   (- 4.584407128516051d0 0)
;;   .1 0))

(defun build-marginal-tree (time-matrix n-number-of-species list-of-descending-coal-times)
  (if ))

1. Make a list of coalescent times in ascending order for each leaf
2. Do we also include the leaf start times?

(defun get-matrix-entry (matrix row column)
  "Retrieve an entry of an array as if it were a matrix (indexing starts at 1)."
  (aref matrix (1- row) (1- column)))

;; ((4.584407128516051d0 (3) (3) (3) (3))
;;  (4.497471958701756d0 (3) (3) (3) (3))
;;  (4.025203326294506d0 (3) (3) NIL (3))
;;  (2.693733411688349d0 (3) (3) NIL NIL))

(defvar tedges)
(defvar test-edge-set)
(setf tedges '(((4.584407128516051d0 (2 1 3) (2 1 3) (1 2 3) (1 2 3)))
              ((4.584407128516051d0 (2 1 3) (2 1 3) (1 2 3) (1 2 3))
               (4.497471958701756d0 (3) (3) (3 2 1) (3))
               (4.025203326294506d0 (1 2 3) (1 2 3) NIL (3 2 1))
               (2.693733411688349d0 (3 2 1) (3 2 1) NIL NIL))))
(setf test-edge-set (second tedges))


(defun replace-entry (matrix row column new-value)
"Replaces a specified matrix element with a new value. Indices start with 1. The function also returns the updated value. Actually uses array data structure."
  (setf
   (apply #'aref matrix (list (1- row) (1- column)))
   new-value))


(defvar *time-matrix*)
(defun construct-time-matrix (n-number-of-species edge-set site-number)
  "INPUT: an edge set (i.e. (second output-edges) or (union (first output-egdes) (second output-edges))) as well as the number of species and a site number. OUTPUT: a matrix of TIMES of coalescences. Let this matrix be M=(m_{ij}). Then 2m_{ij} - d_i - d_j is the coalescence time of i and j, where d_i is the START time of loci i"
  (loop
    initially
       (setf *time-matrix*
             (make-array (list n-number-of-species n-number-of-species)
                                       :initial-element -1))
    for i from 1 to (1- n-number-of-species)
    do (loop for j from (1+ i) to
             n-number-of-species
             do (replace-entry *time-matrix* i j
                               (get-tmrca i j edge-set site-number))))
  *time-matrix*)


(defun get-tmrca (species1 species2 edge-set site-number)
  (loop for x in edge-set
        if (common-coordinate species1 species2 x site-number)
          minimizing (first x)))

   
  ;; know the leaf times and the root time
  ;; know the number of leaves
  ;; need a tree with the associated parameters :dist-from-parent and :mutation-rate
  ;; 

;; Intervals
;;
;; Here we use intervals of the form (m . n) = {m,m+1,...,n}, where m=<n.
;;
;; data format will be a list of disjoint intervals which are ordered increasing
;; by magnitude. Call this a diset (for "disjoint interval set"). For example:
;; '((1 . 5) (6 . 10) (12 . 29)). Nonexample: '((1. 5) (12 . 29) (6 . 10))
;;
;; Useful facts: the union of overlapping intervals is an interval, and that two
;; intervals (a . b) and (c . d) overlap iff a=<d and c=<b

(defun intervalp (x)
  "Test if x is an integer interval."
  (and (consp x)
       (integerp (car x))
       (integerp (cdr x))
       (<= (car x) (cdr x))))

(defun disetp (x)
  "Test if x is a list of integer intervals."
  (and (every 'intervalp x)
       (disjointp x)))

(defun disjointp (x)
  "Test if the intervals in an iset are pairwise disjoint"
  (if (null (rest x))
      t
      (if (< (upper (first x)) (lower (second x)))
          (disjointp (rest x))
          nil)))

(defun overlapp (I J)
  "Test whether two integer intervals overlap."
  (and (<= (car I) (cdr J))
       (<= (car J) (cdr I))))

(defun make-interval (a b)
  "Make an interval with endpoints a and b."
  (declare (type integer a b))
  (if (<= a b)
      (cons a b)
      (cons b a)))

(defun interval-union (I J)
  "unite two intervals"
  (if (overlapp I J)
      (make-interval (min (car I)
                          (car J))
                     (max (cdr I)
                          (cdr J)))
      (list I J)))

(defun interval-overlapping-union (list-of-intervals)
  "unite any number of overlapping intervals"
  (make-interval (loop for x in list-of-intervals minimizing (car x))
                 (loop for x in list-of-intervals maximizing (cdr x))))

(defun find-overlapping-intervals (I D)
  (loop for x in D if (overlapp x I) collecting x))

(defparameter *diset* '((1 . 4) (6 . 10) (12 . 40)))

(defun add-interval-to-diset (I D)
  "Input: an interval I and a diset D. Output: a diset D' consisting of the
union of I and D. Example usage: (add-interval-to-diset '(1 . 111) '((1 . 2) (3
. 11) (12 . 15)))"
  (loop for x in D
        if (overlapp x I)
          collect x into overlapping
        else
          collect x into not-overlapping
        finally
           (return (cons
                    (interval-overlapping-union
                     (cons I overlapping))
                    not-overlapping))))

(defun upper (I)
  (cdr I))
(defun lower (I)
  (car I))


(defun split-diset (breakpoint d)
  "Return a list (left-diset right-diset) consisting of those points in the
diset d which are strictly to the left of the breakpoint and those which
are (not)strictly to the right. Example useage: (split-diset 4 '((1 . 3) (5 .
7)))"
  (loop for x in d
        if (< (upper x) breakpoint)
          collect x into left-diset
        else
          if (<= (lower x) breakpoint)
            collect (make-interval (lower x)
                                   (1- breakpoint))
            into left-diset and
            collect (make-interval breakpoint
                                    (upper x))
            into right-diset
        else
          collect x into right-diset
        finally (return (list left-diset right-diset))))
        
;; this works pretty good, but we need to deal with uniting e.g. the last two
;; intervals in the following: '(((1 . 1) (3 . 3)) ((4 . 5) (6 . 7)))

;; there are two cases: either the breakpoint is contained in one of the
;; intervals or it is not. If it is not, we can just sort the intervals using
;; loop.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
