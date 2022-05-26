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

(defun make-leaf (leaf-label distance-from-parent mutation-rate recombination-rate)
  "Create a leaf. Leaf-label should be a string. The other input parameters
should be numbers."
  `(,leaf-label
    ((:dist-from-parent . ,distance-from-parent)
     (:mutation-rate . ,mutation-rate)
     (:recomb-rate . ,recombination-rate))))

(defun make-node (node-label distance-from-parent mutation-rate
                  recombination-rate left-subtree right-subtree)
  "Create a node. Node-label should be a string, distance-from-parent,
mutation-rate, and recomb-rate should be numbers, and the subtrees subtrees
should be trees (i.e., either nodes or leafs, but not nil)."
  `((,node-label
     ((:dist-from-parent . ,distance-from-parent)
      (:mutation-rate . ,mutation-rate)
      (:recomb-rate . ,recombination-rate)))
    ,left-subtree
    ,right-subtree))

;;______________________________________________________________________________
;;
;; Part 2. Recognizers
;;______________________________________________________________________________

(defun leafp (tree)
  "Return t if tree is a leaf and nil otherwise. Does not distinguish between
leaves with/without vertex times."
  (and (listp tree)
       (= 2 (length tree))
       (stringp (first tree))
       (listp (second tree))))
  
(defun nodep (tree)
  "Return t if tree is a node and nil otherwise. Does not distinguish between
  leaves with/without vertex times."
  (and (listp tree)
       (= 3 (length tree))
       (or (nodep (second tree))
           (leafp (second tree)))
       (or (nodep (third tree))
           (leafp (third tree)))))

;;______________________________________________________________________________
;;
;; Part 2. Selectors
;;______________________________________________________________________________

(defmacro get-parameter (parameter tree)
  "Retrieve specific parameters from a tree node. Example usage: (get-parameter
:mutation-rate *nv1*)"
  `(cdr (assoc ,parameter (get-parameter-alist ,tree))))
(defun get-parameter-alist (tree)
  "Retrieve the parameter alist for the vertex at the top of the input tree.
Works for both leaves and nodes."
  (if (leafp tree)
      (second tree)
      (second (first tree))))
(defun left-subtree (tree)
  "Return the left subtree of tree. If tree is a leaf, return nil."
  (if (leafp tree) nil (second tree)))
(defun right-subtree (tree)
  "Return the right subtree of tree. If tree is a leaf, return nil."
  (if (leafp tree) nil (third tree)))
(defun get-vertex-name (tree)
  "Return the name (i.e. label) of vertex. Works for both leaves and internal
nodes. For leaves, it returns the leaf label. For nodes, it returns the name of
the node."
  (first (if (leafp tree) tree (first tree))))

;; (defun parent-test (tree child-name)
;;   "Test if tree is the parent of the subtree with top vertex labeled child-name."
;;   (or (equal (get-vertex-name (left-subtree tree))
;;              child-name)
;;       (equal (get-vertex-name (right-subtree tree))
;;              child-name)))
      
;; (defun get-parent (tree child-name)
;;   "Returns the parent node of the vertex with label child-name, which is a
;; string. Usage: (get-parent *n3* ``ab'')."
;;   (cond ((leafp tree) nil)
;;         ((parent-test tree child-name) tree)
;;         (t (append (get-parent (left-subtree tree) child-name)
;;                    (get-parent (right-subtree tree) child-name)))))

(defun get-parent (tree child-name)
  "Returns the parent node of the vertex with label child-name, which is a
string. Usage: (get-parent *n3* ``ab''). To understand this recursive function,
one needs to understand two properties of cond: first, cond terminates when the
first condition (evaluated in order) evaluates to a non-nil result. Second, if
test1 is nonnil, then a cond condition of the form ((test1)) returns the value
of test1."
  (cond ((leafp tree) nil)
        ((or (equal (get-vertex-name (left-subtree tree)) child-name)
             (equal (get-vertex-name (right-subtree tree)) child-name))
         tree)
        ((get-parent (left-subtree tree) child-name))
        ((get-parent (right-subtree tree) child-name))))

;;______________________________________________________________________________
;;
;; Part 3. Add vertex times to the trees
;;______________________________________________________________________________

;; The next function counts the number of leaves on a binary tree.
(defun count-number-of-leaves (tree)
  "Count the number of leaves on a binary tree."
  (if (leafp tree)
      1
      (+ (count-number-of-leaves (left-subtree tree))
         (count-number-of-leaves (right-subtree tree)))))

;; The next function computes the total height of a tree.
(defun get-tree-height (tree &optional (x nil))
  "Return the total height of the tree, i.e. from the top of the root edge to
  the maximally distant leaf. Note: the optional variable x is an auxilliary
  variable used to add up the total height of the tree over recursive calls; it
  is used for the recursion and should not be specified by the user.
  Usage: (get-tree-height tree)"
  (if (null x)
      (get-tree-height tree 0)
      (if (leafp tree)
          (+ x (get-parameter :dist-from-parent tree))
          (max (get-tree-height (right-subtree tree)
                (+ x (get-parameter :dist-from-parent tree)))
               (get-tree-height (left-subtree tree)
                (+ x (get-parameter :dist-from-parent tree)))))))
                               
;; We need the following auxilliary function.
(defun add-age-parameters-to-leaf (tree parent-age)
  "Cons (add) the population start and end times to the parameter alist of a leaf object. Note that the population-start-time is the age of the vertex, and the population-end-time is the age of its parent."
  (list (get-vertex-name tree)
        (acons :population-end-time parent-age
               (acons :population-start-time
                      (- parent-age (get-parameter :dist-from-parent tree))
                      (get-parameter-alist tree)))))

;; The next function computes vertex times for a tree and outputs a new tree
;; containing vertex times in the vertex parameter lists.
(defun add-age-parameters-to-tree (tree
                                   &optional (parent-age
                                              (get-tree-height tree)))
  "Input: a tree constructed with the functions 'make-node' and 'make-leaf'.
Output: the same tree, but with labelled vertex times, in the form of two
parameters: population-start-time and population-end-time. Here we regard the
edge extending up from the vertex to be a 'population.' Hence, the
population-start-time is the age of the vertex, and the population-end-time is
the age of its parent. The optional variable is used for recursion and should
not be entered by the user. Usage: (add-age-parameters-to-tree tree)"
  (if (leafp tree)
      (add-age-parameters-to-leaf tree parent-age)
      (list (add-age-parameters-to-leaf tree parent-age)
            (add-age-parameters-to-tree
             (left-subtree tree)
             (- parent-age (get-parameter :dist-from-parent tree)))
            (add-age-parameters-to-tree
             (right-subtree tree)
             (- parent-age (get-parameter :dist-from-parent tree))))))

;; The next function adds a root label to the parameter list
(defun add-root-label (tree)
  "Add a root indicator to the parameter-alist of a tree"
  (when (nodep tree)
    (list
     (list
      (get-vertex-name tree)
      (acons :rootp t (get-parameter-alist tree)))
     (left-subtree tree)
     (right-subtree tree))))


;; The next function adds a root label to the parameter list
(defun add-root-label (tree)
  (list
   (list
    (get-vertex-name tree)
    (acons :rootp t (get-parameter-alist tree)))
   (left-subtree tree)
   (right-subtree tree)))


;;______________________________________________________________________________
;;
;; Part 4. Tree examples for testing
;;______________________________________________________________________________
;; From here on we will assume that all trees have vertex times, i.e. were
;; processed with the function 'add-age-parameters-to-tree'

;; Example 1. A leaf
(defparameter *l1*
  (make-leaf "A" .1 1 4))
(defparameter *lv1*
  (add-age-parameters-to-tree *l1*))

;; Example 2. A node (unbalanced quartet):
(defparameter *n1*
  (make-node "root" 0 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv1* (add-age-parameters-to-tree *n1*))

;; Example 3. A node (unbalanced quartet):
(defparameter *n2*
  (make-node "root" 999 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv2* (add-age-parameters-to-tree *n2*))

;; Example 4. A node (balanced quartet):
(defparameter *n3*
  (make-node "root" 999 .1 .2
             (make-node "ab" 2 .1 .2
                        (make-leaf "a" 2 .1 .2)
                        (make-leaf "b" 2 .1 .2))
             (make-node "cd" 1 .1 .2
                        (make-leaf "c" 3 .1 .2)
                        (make-leaf "d" 2 .1 .2))))
(defparameter *nv3* (add-age-parameters-to-tree *n3*))

;; Example 5. Three taxa example
(defparameter *3-taxa-example* '(("abc" (4 999 1 0))
                                 (("ab" (2 2 0.1 0.2))
                                  ("a" (0 2 0.1 0.2))
                                  ("b" (0 2 0.1 0.2)))
                                 ("c" (0 4 1 0))))

(defparameter *t* '(("abc" (0 1 0))
                    (("ab" (.9 0.1 0.2))
                     ("a" (.5 0.1 0.2))
                     ("b" (.4 0.1 0.2)))
                    ("c" (.8 1 0))))

(defparameter *bigbad*
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
(defparameter *bigbad* (add-age-parameters-to-tree *bigbad*))
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
        

(defun implement-substitution (current-nucleotide &optional vertex-name)
  "Implements a substitution to a *different* nucleotide. Specifically, return a
random number from the set {0,1,2,3}\{current-nucleotide}. It does this by
adding 1d3 to the current nucleotide and then reducing the result modulo 4."
  (declare (integer current-nucleotide))
  (let ((new-nucleotide (mod (+ current-nucleotide 1 (random 3)) 4)))
    (format t "~%Substitution from ~a to ~a in population ~a~%"
            current-nucleotide
            new-nucleotide
            vertex-name)
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
  (let* ((vertex-name (get-vertex-name tree))
         (edge-length (get-parameter :dist-from-parent tree))
         (mutation-rate (get-parameter :mutation-rate tree))
         (substitution-probability
           (* .75 (- 1 (exp (* (/ -4 3) mutation-rate edge-length)))))
         (new-nucleotide (if (< (random 1d0) substitution-probability)
                             (implement-substitution parent-nucleotide vertex-name)
                             (progn
                               (format t "~%No substitution in population ~a~%"
                                       vertex-name)
                               parent-nucleotide))))
    (progn
      (if (nodep tree)
          (unless (equal vertex-name "root")
            (format t "ancestor ~a: ~a ~%" vertex-name new-nucleotide))
          (format t "leaf ~a: ~a ~%" vertex-name new-nucleotide))
      (list vertex-name new-nucleotide))))

;; The next two functions are the workhorses of our sequence simulator. They use
;; the above functions to recursively implement the JC69 process on the input
;; tree.

(defun evolve-down-tree (tree)
  "Recursive implmentation of JC69 on the tree. Starts with random nucleotide at
the tip of the edge above the root. To make this start exactly at the root (by
which I mean the mrca of the samples) you should set the 'dist-from-parent' in
the root vertex to zero."
  (let ((root-state (draw-random-nucleotide)))
    (progn (format t "Root state: ~a ~%" root-state)
           (evolve-down-tree-aux tree root-state nil))))

(defun evolve-down-tree-aux (tree parent-nucleotide x)
  "x is the accumulator variable use to construct the output list"
  (if (leafp tree)
      (cons (evolve-down-edge tree parent-nucleotide) x)
      (let ((new-parent-nucleotide
             (second (evolve-down-edge tree parent-nucleotide))))
        (evolve-down-tree-aux (left-subtree tree)
                              new-parent-nucleotide
                              (evolve-down-tree-aux (right-subtree tree)
                                                    new-parent-nucleotide
                                                    x)))))

;; It works! Example usage:
;; (evolve-down-tree *n1*)

;; The following function allows us to avoid the strange backticks like those
;; used in our original simulator, simulate-three-species

(defun make-leaf-sample (i n-total-leaf-number k-sequence-length leaf)
  "makes the i-th initial sample for a phylogenetic tree with a number of taxa
equal to n-total-leaf-number, and when sequences are length k-sequence-length"
  (progn
    (format t "~%Sample ~a created in population ~a at time ~a~%"
            i
            (get-vertex-name leaf)
            (get-parameter :population-start-time leaf))
    (list
     (list
      (cons
       (get-parameter :population-start-time leaf)
       (append (make-list (- i 1))
               (cons (interval 1 k-sequence-length)
                     (make-list (- n-total-leaf-number i)))))))))

(defun merge-output-edge-sets (child-1 child-2)
  "Combine the output edge sets child-1=(P_1,Q_1) and child-2=(P_2,Q_2) from two
daughter populations for entry into their parent population. Output a pair
of (P,Q) of active and inactive edges in the appropriate format for use in
arg-builder. Note that if one of the edge-sets is nil, this will just output the
other edge set."
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
(load "/home/ultralisk/MSCR-simulator/scripts/simulator.lisp")


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
  (progn
    (defparameter *leaf-number-tracker* 0) ; can I change this to a let form?
    (mscr (count-number-of-leaves species-tree)
          k-sequence-length
          species-tree
          nil)))

(defparameter *leaf-number-tracker* 0)
(defun mscr (n-total-leaf-number k-sequence-length tree edge-sets)
  "Input comments: the variable 'n-total-leaf-number' is number of leaves on the
initial (i.e. full) tree, 'k-sequence-length' is just the length of the
sequences in base pairs. The input 'leaf-number-tracker' tracks the leaf number:
its initial value is 1? and it increments leaf-number-tracker by 1 each time a
right subtree is recursively evaluated. The initial value of edge-sets should be
nil, I think."
  (build-single-population-arg
   (get-parameter :recomb-rate tree)           ; ρ
   (get-parameter :population-start-time tree) ; t_0
   (get-parameter :population-end-time tree)   ; t_end
   (get-vertex-name tree)                      ; population-name
   k-sequence-length
   (if (leafp tree)
       (make-leaf-sample (incf *leaf-number-tracker*)
                         n-total-leaf-number
                         k-sequence-length
                         tree)
       (merge-output-edge-sets (mscr n-total-leaf-number
                                     k-sequence-length
                                     (left-subtree tree)
                                     edge-sets) 
                               (mscr n-total-leaf-number
                                     k-sequence-length
                                     (right-subtree tree)
                                     edge-sets)))
   (equal "root" (get-vertex-name *nv1*)))) ; need to test this -- is it
                                            ; terminating correctly?



(defun construct-marginal-gene-tree (site-number species-tree output-edges)
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
