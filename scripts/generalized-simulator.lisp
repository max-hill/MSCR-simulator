;;;; generalized-simulator.lisp --- Implement the MSCR-JC(k) process on
;;;;                                arbitrary binary trees. Currently under
;;;;                                construction.
;;
;; Author: max hill 
;; (Last updated 2022-06-05)

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
;; 4. (low priority) sensible naming of variables and things
;;    I should rename 'edge' to 'lineage'. Also consider using p-lists:
;;    '(:age 0.1 :ancestral-sites-from-taxa-1 ((1 . 4) (6 . 10) (12 . 40))
;;    :ancestral-sites-from-taxa-2 ((1 . 7) (9 . 11))
;;    :ancestral-sites-from-taxa-3 NIL)
;;
;; Sources: The following resource was extremely helpful:
;; https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html
;;

;; We start by defining the following global parameters. If *mscr-verbose-mode*
;; is non-nil, then the MSCR simulator will print output about what it is doing.
;; This is for testing purposes. 
(defparameter *mscr-verbose-mode* nil)
(defparameter *mutation-rate* 1)

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



(defun make-leaf
    (&key
       (population-name "leaf")
       (dist-from-parent nil)
       (mutation-rate *mutation-rate*)
       (recomb-rate 0)
       (parent-label "parent")
       (population-start-time nil)
       (population-end-time nil)
       (numeric-label nil))
  "Create a leaf population (i.e. a population corresponding to a leaf edge of
the species tree. Consists of a list containing a single element: an alist of
parameters."
  (progn
    (cond ((and (not population-end-time) dist-from-parent population-start-time)
           (setf population-end-time (+ population-start-time dist-from-parent)))
          ((and (not dist-from-parent) population-start-time population-end-time)
           (setf dist-from-parent (- population-end-time population-start-time)))
          ((and (not population-start-time) population-end-time dist-from-parent)
           (setf population-start-time (- population-end-time dist-from-parent))))
    (list
     (list :population-name population-name
           :dist-from-parent dist-from-parent
           :mutation-rate mutation-rate
           :recomb-rate recomb-rate
           :parent-label parent-label
           :population-start-time population-start-time
           :population-end-time population-end-time
           :numeric-label numeric-label))))

(defun make-node
    (&key
       (population-name "node")
       (left-subtree nil)
       (right-subtree nil)
       (dist-from-parent nil)
       (mutation-rate *mutation-rate*)
       (recomb-rate 0)
       (parent-label "parent")
       (population-start-time nil)
       (population-end-time nil)
       (numeric-label nil))
  "Create a node population (i.e. a population corresponding to an internal node
and edge of the species tree. Consists of a list containing three elements: an
alist of parameters, a left subtree, and a right subtree."
  (cons
   (first
    (make-leaf :population-name population-name
              :dist-from-parent dist-from-parent
              :mutation-rate mutation-rate
              :recomb-rate recomb-rate
              :parent-label parent-label
              :population-start-time population-start-time
              :population-end-time population-end-time
              :numeric-label numeric-label))
   (list left-subtree
         right-subtree)))

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
;; Part 3. Selectors
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
  "Return the name (i.e. label) of the population. Works for both leaves and
internal nodes."
  (get-parameter :population-name tree))

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
        
;;______________________________________________________________________________
;;
;; Part 4. Augmenting tree Parameters
;;______________________________________________________________________________
;; This section defines a function which computes additional information about
;; the tree and adds it to the parameter list.

;; First we define some counting functions to obtain the heighet as well as the
;; number of leaves on the binary tree.
(defun count-number-of-leaves (tree)
  "Count the number of leaves on a binary tree."
  (if (leafp tree)
      1
      (+ (count-number-of-leaves (left-subtree tree))
         (count-number-of-leaves (right-subtree tree)))))

(defun get-tree-height-aux (tree x)
  "Auxilliary function for get-tree-height. The variable x is a dummy variable
  for adding up the total height of the tree."
  (let ((dist-from-parent (get-parameter :dist-from-parent tree)))
    (if (leafp tree)
        (+ x dist-from-parent)
        (max
         (get-tree-height-aux (right-subtree tree) (+ x dist-from-parent))
         (get-tree-height-aux (left-subtree tree) (+ x dist-from-parent))))))

(defun get-tree-height (tree)
  "Return the total height of the tree, i.e. from the top of the root edge to
  the maximally distant leaf."
  (get-tree-height-aux tree 0))

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
         (current-population-name (get-parameter :population-name tree)))
    (set-population-parameter :parent-label parent-label tree)
    (set-population-parameter :population-start-time current-population-start-time tree)
    (set-population-parameter :population-end-time parent-age tree)
    (unless (leafp tree)
      (add-ages-and-parents (left-subtree tree)
                            current-population-start-time
                            current-population-name)
      (add-ages-and-parents (right-subtree tree)
                            current-population-start-time
                            current-population-name))))

;; The next two functions augment the input tree parameter lists with a
;; structured population labeling system.

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

;; Finally, the key function of this section is the following.

(defun augment-tree-parameters (tree)
  "The input should be a tree constructed with the functions 'make-node' and
'make-leaf'. This function is idempotent."
  (progn (add-ages-and-parents tree)
         (add-numeric-labels tree)))

;;______________________________________________________________________________
;;
;; Part 5. Tree examples for testing
;;______________________________________________________________________________
;; From here on we will assume that all trees have vertex times, i.e. were
;; processed with the function 'add-age-parameters-to-tree'

;; Example 1. A leaf
(defparameter a* (make-leaf :population-name "A" :dist-from-parent 1.7 :recomb-rate 1))
(augment-tree-parameters a*)

;; Example tree: an unbalanced quartet:

(defparameter t*
(make-node :population-name "ABCD"
           :dist-from-parent 999
           :recomb-rate 0
           :left-subtree (make-node
                           :population-name "ABC"
                           :dist-from-parent 1
                           :recomb-rate 1.2
                           :left-subtree (make-node
                                           :population-name "AB"
                                           :dist-from-parent 1.3
                                           :recomb-rate 3
                                           :left-subtree (make-leaf
                                                          :population-name "A"
                                                          :dist-from-parent 1.7
                                                          :recomb-rate 1)
                                           :right-subtree (make-leaf
                                                           :population-name "B"
                                                           :dist-from-parent 1.5
                                                           :recomb-rate 0))
                           :right-subtree (make-leaf
                                           :population-name "C"
                                           :dist-from-parent .6
                                           :recomb-rate 1.02))
                           :right-subtree (make-leaf
                                           :population-name "D"
                                           :dist-from-parent 4
                                           :recomb-rate 1.1)))

(augment-tree-parameters t*)

;; Example tree: a balanced quartet:
(defparameter s*
  (make-node :population-name "ABCD"
             :dist-from-parent 999
             :recomb-rate 0
             :left-subtree (make-node
                            :population-name "AB"
                            :dist-from-parent 1
                            :recomb-rate 1.2
                            :left-subtree (make-leaf
                                           :population-name "A"
                                           :dist-from-parent 1.3
                                           :recomb-rate 3)
                            :right-subtree (make-leaf
                                            :population-name "B"
                                            :dist-from-parent 1.5
                                            :recomb-rate 0))
             :right-subtree (make-node
                             :population-name "CD"
                             :dist-from-parent 1.7
                             :left-subtree (make-leaf
                                            :population-name "C"
                                            :dist-from-parent 1
                                            :recomb-rate 1.02)
                             :right-subtree (make-leaf
                                             :population-name "D"
                                             :dist-from-parent 1
                                             :recomb-rate 1.1))))
(augment-tree-parameters s*)

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
    (when *mscr-verbose-mode*
      (format t "~%Substitution from ~a to ~a in population ~a~%"
              current-nucleotide
              new-nucleotide
              population-name))
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
         (substitution-probability
           (* .75 (- 1 (exp (* (/ -4 3)
                               (get-parameter :dist-from-parent tree)
                               (get-parameter :mutation-rate tree))))))
         (new-nucleotide
           (if (< (random 1.0) substitution-probability)
               (implement-substitution parent-nucleotide population-name)
               (progn (when *mscr-verbose-mode*
                        (format t "~%No substitution in population ~a~%" population-name))
                      parent-nucleotide))))
      (when *mscr-verbose-mode*
        (if (leafp tree)
            (format t "leaf ~a: ~a ~%" population-name new-nucleotide)
            (unless (rootp tree)
              (format t "ancestor ~a: ~a ~%" population-name new-nucleotide))))
      (list population-name new-nucleotide)))

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
    (progn (when *mscr-verbose-mode* (format t "Root state: ~a ~%" root-state))
           (evolve-down-tree-aux tree root-state nil))))

;; It works! Example usage:
;; (implement-jukes-cantor-process *n1)

;;______________________________________________________________________________
;;
;; Part 7. Integer Intervals
;;______________________________________________________________________________
;;
;; In this section, we implement an interval system which will be used by the
;; MSCR simulator. In paraticular, we consider sets of integer intervals, where
;; integer intervals take the form (m . n) = {m,m+1,...,n}, where m=<n. We say
;; that two integer intervals (a . b) and (c . d) are separated iff either b+1<d
;; or d+1<a. Obviously, separatedness is strictly stronger than disjointedness.
;; For example, the intervals (1 . 3) and (4 . 8) are disjoint but not
;; separated.
;;
;; If two integer intervals are not separated, we say that they are overlapping.
;; While this terminology is not consistent with regular use of the term
;; 'overlapping', it is useful to give a name to the property.
;;
;; In particular, we will consider lists of separated intervals of integers
;; which are ordered increasing by magnitude. We call such lists osisets (for
;; "ordered separated interval sets". Examples oisets: '((1 . 4) (6 . 10) (12 .
;; 29)). Nonexamples: '((1 . 5) (6 . 10) (12 . 29)) and '((1. 5) (12 . 29) (6 .
;; 10)).
;;
;; Useful facts:
;;   1. The union of overlapping intervals is an interval,
;;   2. Intervals intervals (a . b) and (c . d) overlap iff a=<d+1 AND c=<b+1

;; Basic integer intervals are constructed using the following function.

(defun make-interval (a b)
  "Make an interval with endpoints a and b."
  (declare (optimize (speed 3)))
  (declare (type fixnum a b))
  (if (<= a b)
      (cons a b)
      (cons b a)))

;; The next several functions are the selectors and recognizers for integer
;; intervals and osisets.

(defun upper (I)
  (cdr I))

(defun lower (I)
  (car I))

(defun separatedp (I J)
  "Test if two integer intervals I and J are separated. Setting safety = 0 here
is okay - there is no risk of overflow since we never add these numbers, and
they are bounded by k, which presumably is less than 2^15."
  (declare (optimize (speed 3) (safety 0)))
  (or (< (the fixnum (1+ (the fixnum (cdr I)))) (the fixnum (car J)))
      (< (the fixnum (1+ (the fixnum (cdr J)))) (the fixnum (car I)))))

(defun overlapp (I J)
  "Test if two integer intervals I and J are overlapping."
  (not (separatedp I J)))

(defun intervalp (x)
  "Test if x is an integer interval. The null set is regarded as an interval."
  (or (null x)
      (and (consp x)
           (integerp (car x))
           (integerp (cdr x))
           (<= (car x) (cdr x)))))

(defun osisetp (x)
  "Test if x is an osiset."
  (if (null (rest x))
      (intervalp (first x))
      (and (every 'intervalp x)
           (< (1+ (upper (first x))) (lower (second x)))
           (osisetp (rest x)))))

(defun number-in-intervalp (x interval)
  "Test whether the number x is in the interval (an integer interval)."
  (and (>= x (lower interval))
       (<= x (upper interval))))

(defun test-membership (x osiset)
  "Test whether the integer x is contained in the osiset."
  (loop for interval in osiset
          thereis (and (>= x (car interval))
                       (<= x (cdr interval)))))

;; The next three functions are used for merging osisets. This will be used for
;; implementing coalescences.

(defun unite-overlapping-intervals (list-of-intervals)
  "Input: a list of overlapping intervals. Output: the union of those intervals.
Does not test whether the input intervals are overlapping."
  (make-interval (loop for x in list-of-intervals minimizing (car x))
                 (loop for x in list-of-intervals maximizing (cdr x))))

(defun add-interval-to-osiset (I osiset)
  "Input: an interval I and a osiset D. Output: a osiset consisting of the
union of I and D. Example usage: (add-interval-to-osiset '(1 . 111) '((1 . 2) (3
. 11) (12 . 15)))"
  (loop for x in osiset
        with overlap-indicator = nil
        with intervals-remaining = (length osiset)
        with interval-added = nil
        do (decf intervals-remaining)
        if (overlapp x I)
          collect x into list-of-overlapping-intervals
          and do (setf overlap-indicator t)
          and if (= 0 intervals-remaining)
                collect (unite-overlapping-intervals
                         (cons I list-of-overlapping-intervals))
                  into new-osiset end
        else
          if overlap-indicator
            do (setf overlap-indicator nil interval-added t)
            and collect (unite-overlapping-intervals
                         (cons I list-of-overlapping-intervals))
                  into new-osiset
            and collect x into new-osiset
        else
          if (and (not interval-added) (< (upper I) (lower x)))
            collect I into new-osiset and do (setf interval-added t) end
          and collect x into new-osiset
          and if (and (> (lower I) (upper x)) (= 0 intervals-remaining))
                collect I into new-osiset end
                finally (return new-osiset)))                  

(defun merge-osisets (osiset1 osiset2)
  "Inpute: two osisets. Output: a new osiset which equals the union of the two
input osisets."
  (cond ((endp osiset1) osiset2) ; if either osiset is empty, just return 
        ((endp osiset2) osiset1) ; the other one.
        (t (loop for interval in osiset1
                 with new-osiset = osiset2
                 do (setf new-osiset (add-interval-to-osiset interval new-osiset))
                 finally
                    (return new-osiset)))))

;; The next function splits osisets. This will be used for implementing
;; recombinations.

(defun split-osiset (breakpoint input-osiset)
  "Return a list (left-osiset right-osiset) consisting of those points in the
osiset d which are strictly to the left of the breakpoint and those which
are (not)strictly to the right. Example useage: (split-osiset 4 '((1 . 3) (5 .
7)))"
  (loop for x in input-osiset
        if (< (upper x) breakpoint) collect x into left-osiset
          else
            if (<= breakpoint (lower x)) collect x into right-osiset
              else
                collect (make-interval (lower x) (1- breakpoint)) into left-osiset
                and collect (make-interval breakpoint (upper x)) into right-osiset
        finally (return (list left-osiset right-osiset))))

;;______________________________________________________________________________
;;
;; Part 8. MSCR Simulator
;;______________________________________________________________________________
;;
;; This section implements the MSCR simulator on a general binary tree.

;; Next we define several general auxillary functions. Most of them were copied
;; over from simulator.lisp without change.

(defun draw-exponential (λ)
  "Return a number drawn according to a rate λ exponential random variable.
Based on code from https://github.com/tpapp/cl-random. If λ=0, return positive
infinity, represented by most-positive-long-float (1.7976931348623157d308)."
  (if (zerop λ)
      most-positive-long-float
      (- (/ (log (- 1 (random 1d0))) λ))))

(defun randomly-choose (x &optional (n 1) (with-replacement nil))
  "Randomly choose n elements from a list x without (or with) replacement. If
the number n is unspecified (or n=1), output a single randomly-chosen element.
If n>1, output a *list* of elements. By default, the selection of multiple
elements is done *without* replacement. Selection *with* replacement can be
specified by setting the optional variable 'with-replacement' to a non-nil
value."
  (cond ((null x) nil)
	((= n 1) (nth (random (length x)) x))
	(t (randomly-choose-several x n with-replacement))))

(defun randomly-choose-several (x &optional (n 1) (with-replacement nil))
  "This is an auxillary function to the function randomly-choose. It outputs a
*list* of n elements chosen randomly without (or with) replacement from the list
x. It can also be run standalone. The only difference in output compared to
randomly-choose occurs when n=1, in which case randomly-choose outputs only the
randomly-chosen element, whereas randomly-choose-several outputs a list of
length one containing the element. I don't know if that will ever be of
interest, but I record it here just in case."
  (cond ((null x) nil)
        ((<= n 0) nil)
	(with-replacement
	    (cons (randomly-choose x)
		  (randomly-choose-several x (1- n) with-replacement)))
	((>= n (length x)) x)
        (t (let ((element (randomly-choose x)))
	     (cons element
		   (randomly-choose-several
		    (remove element x :count 1)
		    (1- n)))))))

(defun remove-elements (elements-to-remove initial-set)
  "Code based on function found at
http://www-users.cselabs.umn.edu/classes/Spring-2018/csci4511/lisp/mapping.html"
  (remove-if
   #'(lambda (element) (member element elements-to-remove :test #'equal))
   initial-set))

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

(defun get-active-lineages (edge-sets)
  "For code readability. Input: (P,Q). Returns P"
  (first edge-sets))

;; Next we define several auxillary functions needed for implementing specific
;; steps in the MSCR process. These have been modified from the versions in
;; simulator.lisp to (1) allow for general tree input rather than just trees
;; with three leaves, and (2) to allow for symbolic interval system, rather than
;; just tracking sites using the function interval (which is too slow).

(defvar *list-of-breakpoints* nil)

(defun make-recombination-parents (time edge breakpoint)
  "Example: (make-recombination-parents .99 `(.1 ,*osiset ((1 . 7) (9 . 11)) nil) 3)"
  (let ((list-of-osisets (rest edge)))
    (loop for osiset in list-of-osisets
          for split-set = (split-osiset breakpoint osiset)
          collect (first split-set) into left-parent
          collect (second split-set) into right-parent
          finally (return (list (cons time left-parent) (cons time right-parent))))))


(defun implement-recombination (time edge-sets number-of-base-pairs
                                &optional (population-name nil))
  "Updates the edge-sets (p,q) appropriately for when a coalescence occurs at
the given time."
  (let* ((recombination-child (randomly-choose (first edge-sets)))
	 (breakpoint (+ 2 (random (- number-of-base-pairs 1)))) ; if x < break --> x goes left
                                                                ; if x >= break -> x goes right
                                                                ; so this choice of +2 is correct.
	 (recombination-parents (make-recombination-parents time recombination-child breakpoint))
	 (new-p (cons (first recombination-parents)
		      (cons (second recombination-parents)
			    (remove-elements
			     (cons recombination-child nil)
			     (first edge-sets))))))
    (progn
      (setf *list-of-breakpoints* (cons breakpoint *list-of-breakpoints*))
      (when *mscr-verbose-mode*
        (format t "~%~%RECOMBINATION in ~a at time ~a" population-name time)
        (format t "~%Breakpoint: ~a" breakpoint)
        (format t "~%One child edge removed: ~a" recombination-child)
        (format t "~%Two parent edges added: ~a" (first recombination-parents))
        (format t "~%                        ~a" (second recombination-parents)))
      (list new-p (second edge-sets)))))

;; For making initial an linages sets (P,Q) at the start of each leaf
;; population, we use the following function.

(defun make-leaf-sample (i n-total-leaf-number k-sequence-length leaf)
  "Make the i-th initial sample for a phylogenetic tree with a number of taxa
equal to n-total-leaf-number, and when sequences are length k-sequence-length."
  (let ((new-lineage
          (list
           (list
            (cons (get-parameter :population-start-time leaf)
                  (append (make-list (- i 1))
                          (cons (list (make-interval 1 k-sequence-length))
                                (make-list (- n-total-leaf-number i)))))))))
    (progn
      (when *mscr-verbose-mode*
        (format t "~%~%LINEAGE ~a CREATED in ~a at time ~a~%One lineage added: ~a"
                i (get-population-name leaf) (get-parameter :population-start-time leaf)
                (first (first new-lineage))))
      new-lineage)))
;; example: (make-leaf-sample 1 4 5 *lv1) remember this gives a pair (P,Q) in
;; which Q=nil and P is a singleton containing only one edge. So the edge is
;; actually (first (first (make-leaf-sample1 1 4 5 *lv1)))

(defun make-coalescent-parent (time coalescing-pair)
  "Creates parent edge of two coalescing edges. The input coalescing-edges is of
the form (x y) where x and y are the edges"
  (let* ((edge1 (first coalescing-pair))
	(edge2 (second coalescing-pair))
        (list1-of-osisets (rest edge1)) ; rename osiset to 'ancestor set'?
        (list2-of-osisets (rest edge2))) ; rename osiset to 'ancestor set'?
    (cons time
          (loop for osiset1 in list1-of-osisets
                for osiset2 in list2-of-osisets
                collecting (merge-osisets osiset1 osiset2) into x
                finally (return x)))))

;; This implements coalescences.

(defun implement-coalescence (time edge-sets &optional (population-name nil))
  "Updates the edge-sets (p,q) appropriately for when a recombination occurs at
the given time."
  (let* ((active-lineages (first edge-sets))
         (inactive-lineages (second edge-sets))
         (coalescing-pair (randomly-choose active-lineages 2))
         (coalescent-parent (make-coalescent-parent time coalescing-pair))
         (new-p (cons coalescent-parent (remove-elements coalescing-pair active-lineages)))
         (new-q (cons coalescent-parent inactive-lineages)))
    (progn
      (when *mscr-verbose-mode*
        (format t "~%~%COALESCENCE in ~a at time ~a" population-name time)
        (format t "~%Two child edges removed: ~a" (first coalescing-pair))
        (format t "~%                         ~a" (second coalescing-pair))
        (format t "~%One parent edge created: ~a" coalescent-parent))
      (list new-p new-q))))

;; The next function is the key auxillary function to the main simulator; it
;; simulates the ARG process in a *single* population. It is analogous to
;; 'arg-builder' in simulator.lisp.

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

;; The next two functions define the top-level of the MSCR simulator.

(defun mscr-aux (n k tree)
  "Auxilliary function for mscr. The variable 'n' is number of leaves on the
initial (i.e. full) tree, 'k' is the length of the sequences in base pairs."
  (if (leafp tree)
      (make-leaf-sample (get-parameter :numeric-label tree) n k tree)
      (let* ((lineages-from-left-daughter (mscr-aux n k (left-subtree tree)))
             (lineages-from-right-daughter (mscr-aux n k (right-subtree tree)))
             (starting-lineages (merge-output-edge-sets lineages-from-left-daughter
                                                        lineages-from-right-daughter)))
        (build-single-population-arg (get-parameter :recomb-rate tree)
                                     (get-parameter :population-start-time tree)
                                     (get-parameter :population-end-time tree)
                                     (get-parameter :population-name tree) ; maybe use the numeric labels?
                                     k
                                     starting-lineages
                                     (rootp tree)))))

(defun mscr (species-tree k-sequence-length)
  (progn
    (defparameter *list-of-breakpoints* nil)
    (mscr-aux (count-number-of-leaves species-tree) k-sequence-length species-tree)))

;; We could simplify the code even further by making k,n, and edge-sets global
;; variables. In particular, if edge-sets were a global variable, we could
;; modify it using setf. I don't think it would be appreciably faster, but the
;; code would probably be shorter.

;;______________________________________________________________________________
;;
;; Part 9. MISC from Part 8: notes, vars, examples, and alternative functions.
;;______________________________________________________________________________
;;

(defparameter *leaf-sample1 (make-leaf-sample 1 4 500 t*))
(defparameter *leaf-sample2 (make-leaf-sample 2 4 500 t*))
(defparameter *leaf-sample3 (make-leaf-sample 3 4 500 t*))
(defparameter *leaf-sample4 (make-leaf-sample 4 4 500 t*))
(defparameter test1* (list (interval 1 30) (interval 1 30)))
(defparameter test2* (list (interval 31 60) (interval 31 60)))

(defparameter *osiset1 '((1 . 4) (6 . 10) (12 . 15) (19 . 20)))
(defparameter *osiset2 '((1 . 3) (9 . 13) (25 . 30)))
(merge-osisets *osiset1 *osiset2)

(defparameter *osiset '((1 . 4) (6 . 10) (12 . 40)))
(defparameter *nonexample '((1 . 5) (6 . 10) (12 . 40)))
(defparameter *osiset '((1 . 2) (4 . 59) (61 . 65) (70 . 80) (90 . 100)))

(defparameter tedges '(((4.584407128516051d0
                           ((2 . 3) (5 . 10))
                           ((4 . 4) (7 . 10))
                           ((1 . 10))))))

;; Old format for input edges: 
;; INPUT EDGE: (.1 ,(interval 1 7)
;;                 ,(interval 3 6)
;;                 ,(interval 2 10))
;; SAMPLE EDGE: (NIL NIL (1 2 3 4 5 6 7) NIL NIL)

;; New format for input edges:
;; INPUT EDGE: (.1 osiset-1 ... osiset-n)

;; examples:
(make-coalescent-parent 999 (list (first (first *leaf-sample1)) (first (first *leaf-sample2))))
(defparameter edge1* (first (first *leaf-sample1)))
(defparameter edge2* '(0.1 nil ((1 . 4) (6 . 10) (12 . 40)) ((1 . 7) (9 . 11)) NIL))
(make-coalescent-parent .99 (list edge1* edge2*))
(merge-output-edge-sets *leaf-sample1 *leaf-sample2)

;; ;; The following verison is somewhat faster. Can use either nconc or append --
;; ;; nconc is faster, I have to check that it doesn't screws anything up.
;; (defun merge-output-edge-sets1 (child-1 child-2)
;;   "Combine the output edge sets child-1=(P_1,Q_1) and child-2=(P_2,Q_2) from two
;; daughter populations for entry into their parent population. Output a pair
;; of (P,Q) of active and inactive edges in the appropriate format for use in
;; arg-builder. Note that if one of the edge-sets is nil, this will just output the
;; other edge set. Recall that the set P consists of those lineages which are
;; active in the sense that they may still undergo coalescent or recombination
;; events. The seq Q contains lineages which have been deactivated as a result of
;; undergoing coalescent events -- but these can't just be thrown out because they
;; contain information about MRCAs."
;;   (list (append (first child-1) (first child-2))
;;         (append (second child-1) (second child-2))))


;; ;; Another version of this file is given here, which does not rely on
;; ;; remove-elements (which might be slow) and does not rely on randomly-choose. I
;; ;; expect this to perform better when n is large, but I could be wrong.
;; (defun implement-coalescence1 (time edge-sets &optional (population-name nil))
;;   (let* ((n (length edge-sets))
;;          (index1 (random n))
;;          (index2 (randomly-choose-excluding n index1))
;;          (active-lineages (first edge-sets))
;;          (coalescing-pair (list (nth index1 active-lineages)
;;                                 (nth index2 active-lineages)))
;;          (coalescent-parent (make-coalescent-parent time coalescing-pair))
;;          (new-p (cons coalescent-parent
;;                       (loop for i from 0 upto (1- n)
;;                             unless (or (= i index1) (= i index2))
;;                               collect (nth i active-lineages))))
;;          (new-q (cons coalescent-parent (second edge-sets))))
;;     ;; (Format T "~%~%COALESCENCE In ~A At Time ~A" Population-Name Time)
;;     ;; (Format T "~%Two Child Edges Removed: ~A" (First Coalescing-Pair))
;;     ;; (Format T "~%                         ~A" (Second Coalescing-Pair))
;;     ;; (Format T "~%One Parent Edge Created: ~A" Coalescent-Parent)
;;     (list new-p new-q)))

;; (defun randomly-choose-excluding (n x)
;;   "Choose a numbers at random from the set {0,...,n-1}\{x}."
;;   (let ((output (random n)))
;;     (if (= output x)
;;         (randomly-choose-excluding n x)
;;         output)))



;;______________________________________________________________________________
;;
;; Part 10. Connecting the MSCR and JC Simulators
;;______________________________________________________________________________



(defun common-coordinate (leaf1 leaf2 edge site)
  "Test whether an edge is ancestral to a specified site sampled from two leaf
samples. In other words, test whether the edge contains that site for both
species."
  (and (test-membership site (nth leaf1 edge))
       (test-membership site (nth leaf2 edge))))
;; COMMENTARY: The inputs 'leaf1' and 'leaf2' take values 1,2,3,
;; corresponding to species A,B,C respectively. To compute all distances, we
;; will need to loop over all pairs (species1,species2) where leaf1 and
;; leaf2 are not equal
;;
;; Example code:
;; * (common-coordinate 1 2 '(3.9 ((1 . 352) (606 . 1000)) ((1 . 352) (606 . 1000)) nil nil) 352)
;; T
;; * (common-coordinate 1 2 '(3.9 ((1 . 352) (606 . 1000)) ((1 . 352) (606 . 1000)) nil nil) 353)
;; NIL

(defun find-descendants (site edge)
  "Return a list of those species which inherit site i from the lineage edge. If
the edge does not have site as an ancestor for any species, return nil."
  (loop for osiset in (rest edge) 
        for species-number from 1
        when (test-membership site osiset)
          collect species-number into descendants
        finally (return descendants)))
;; Example: (find-descendants 23 '(3.9 ((1 . 352) (606 . 1000)) ((1 . 352) (606 . 1000)) nil nil))

;; Consider the following example data:

(defparameter *mscr-output* '(((5.788608309440108d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000))))
                              ((5.788608309440108d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000)))
                               (4.885826667890404d0 ((353 . 605)) ((353 . 605)) ((1 . 1000)) NIL)
                               (4.266150707283983d0 ((1 . 352) (606 . 1000)) ((1 . 352) (606 . 1000)) NIL ((1 . 1000)))
                               (2.949941650510326d0 ((353 . 1000)) ((353 . 1000)) NIL NIL)
                               (1.9481898222457033d0 ((1 . 1000)) ((1 . 1000)) NIL NIL)
                               (3.8150415999214697d0 ((1 . 352)) ((1 . 352)) NIL NIL)
                               (3.9228970527341116d0 ((1 . 352) (606 . 1000)) ((1 . 352) (606 . 1000)) NIL NIL))))
;; which has the following breakpoints
(defparameter *recombination-breakpoints* '(1 352 606))
(defparameter *example-site* 1)

;; We process the data with the following function, which reorders the relevant
;; edges so that the coalescent times are ascending.
(defun process-mscr-output (mscr-output)
  "Return the coalescent lineages, ordered ascending by time. Nondestructive.
Remove copy-list to make it destructive (but faster)."
  (sort (copy-list (second mscr-output)) #'< :key #'car ))

(defparameter *sorted-edges* (process-mscr-output *mscr-output*))

;; The next variable is a list of labels of those lineages which have coalesced.
;; Elements take the form e.g. (1 2 3), which represents the lineage ancestral
;; to leaves 1, 2, and 3. When two lineages coalesce, we merge their
;; corresponding labels by taking a union. When a node coalesces with a leaf, we
;; cons the leaf number on. When two nodes coalesce, we append their labels.
(defparameter *previously-coalesced* '((1 2)))


;; The next variable is an alist of tree components. Elements take the 
(defparameter *tree-builder* nil)

;; PLAN: We will fix a site i, and then loop through *sorted-edges*,
;; constructing a node of the tree for each edge containing site i. Constructed
;; nodes will be stored in *tree-builder*.

;; the following code is for constructing a node when given an input edge and a site i
(progn (defparameter *previously-coalesced* nil)
       (defparameter *tree-builder* nil)
       (loop for edge in *sorted-edges*
             do (progn
                  (format t "~%~%EDGE: ~a" edge)
                  (make-node-from-edge 1 edge))))

(defun already-coalesced-p (r s)
  "Return t if both lineages have already coalesced at a more recent time, which
can be determined by the size of the lists r and s."
  (and (= (length s) 1)
       (= (length r) 0)))

(progn (defparameter *previously-coalesced* nil)
       (defparameter *tree-builder* nil)
       (make-node-from-edge 500 (nth 3 *sorted-edges*))
       (print *previously-coalesced*))

(defun make-node-from-edge (site edge)
  "Analyzes a single edges -- using *previously-coalesced*, determine which lineages have coalesced (wrt a specific site). Update *previously-coalesced* accordingly. And then (not yet implemented) add the corresponding tree node to *tree-builder*."
  (let* ((all-descendants-of-lineage (find-descendants site edge))
         (sr (loop for x in *previously-coalesced*
                with r = all-descendants-of-lineage ;; r = list of remaining lineages
                until (<= (length r) 1)
                when (sort (copy-list (intersection x r)) #'<)
                collect it into s  ;; s = list of coalescing nodes (may have 0, 1 or 2 elements)
                and do (setf r (set-difference r (intersection r x)))
                finally (return (list s r))))
         (s (first sr))
         (r (second sr)))
    (progn (format t "~%sr: ~a~%s: ~a~%r: ~a" sr s r)
           (unless (already-coalesced-p r s)
             (let* ((new-label-to-add (sort (copy-list (union (union (first s) (second s)) r)) #'<)))
               (cond
                 ;; Case 1: both lineages are leaves (base case)
                 ((and (= (length s) 0) (= (length r) 2))
                  (progn
                    (setf *previously-coalesced* (cons new-label-to-add *previously-coalesced*))
                    (format t "~%Both coalescing lineages are leaves.~%New label: ~a" new-label-to-add)
                    ;; Update the tree by adding the node to the alist
                    (setf *tree-builder* (cons
                                          (list new-label-to-add ;this will be the key to the node in the alist
                                                (make-node :population-name new-label-to-add
                                                           :population-start-time (car edge)
                                                           :left-subtree (make-leaf :population-name (first r)
                                                                                    :population-end-time (car edge))
                                                           :right-subtree (make-leaf :population-name (second r)
                                                                                     :population-end-time (car edge))))
                                          *tree-builder*))
                    (format t "~%Updated tree-builder: ~a" *tree-builder*)))
                 ;; Case 2: one coalescing lineage is a node and the other is a leaf
                 ((and (= (length s) 1) (= (length r) 1))
                  (progn
                    (setf *previously-coalesced* (cons new-label-to-add *previously-coalesced*))
                    (format t "~%One coalescing lineage is a leaf; the other is a node.~%New label: ~a~%Label to get node is: ~a"
                            new-label-to-add (first s))
                    (setf *tree-builder* (cons
                                          (list
                                           new-label-to-add
                                           (make-node :population-name new-label-to-add
                                                      :population-start-time (car edge)
                                                      :left-subtree (make-leaf :population-name (car r)
                                                                               :population-end-time (car edge))
                                                      :right-subtree (second (assoc (first s) *tree-builder* :test #'equal))))
                                          *tree-builder*))))
                 ;; Case 3: both coalescing lineages are nodes
                 ((= (length s) 2)
                  (progn
                    (setf *previously-coalesced* (cons new-label-to-add *previously-coalesced*))
                    (format t "~%Both coalescing lineages are nodes.~%New label: ~a~%Labels to get nodes are: ~a ~%~a"
                            new-label-to-add (first s) (second s))
                    (setf *tree-builder* (cons
                                          (list
                                           new-label-to-add
                                           (make-node :population-name new-label-to-add
                                                      :population-start-time (car edge)
                                                      :left-subtree (second (assoc (first s) *tree-builder* :test #'equal))
                                                      :right-subtree (second (assoc (second s) *tree-builder* :test #'equal))))
                                          *tree-builder*))))
		 (t (format t "something unexpected happened"))))))))
                                   

(defun add-times (tree x)
  (cond ((rootp tree)
         (set-population-parameter :dist-from-parent 0 tree))
        ((leafp tree)
         (set-population-parameter :dist-from-parent x tree))
        ((nodep tree)
         (progn (set-population-parameter :dist-from-parent x tree)
                (add-times (left-subtree (get-parameter :populationtree 
                
         
      


((:POPULATION-NAME (1 2 3 4) :DIST-FROM-PARENT NIL :MUTATION-RATE 1
  :RECOMB-RATE 0 :PARENT-LABEL nil :POPULATION-START-TIME
                   5.788608309440108d0 :POPULATION-END-TIME NIL :NUMERIC-LABEL NIL)
 ((:POPULATION-NAME 3 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
   :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME
                    5.788608309440108d0 :NUMERIC-LABEL NIL))
 ((:POPULATION-NAME (1 2 4) :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE
                    0 :PARENT-LABEL "parent" :POPULATION-START-TIME 4.266150707283983d0
   :POPULATION-END-TIME NIL :NUMERIC-LABEL NIL)
  ((:POPULATION-NAME 4 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
    :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME
                     4.266150707283983d0 :NUMERIC-LABEL NIL))
  ((:POPULATION-NAME (1 2) :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE
                     0 :PARENT-LABEL "parent" :POPULATION-START-TIME 1.9481898222457033d0
    :POPULATION-END-TIME NIL :NUMERIC-LABEL NIL)
   ((:POPULATION-NAME 1 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
     :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME
                      1.9481898222457033d0 :NUMERIC-LABEL NIL))
   ((:POPULATION-NAME 2 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
     :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME
                      1.9481898222457033d0 :NUMERIC-LABEL NIL)))))

;; ;; ;; Old program

;; (defun make-node-from-edge (site edge)
;;   "Analyzes a single edges -- using *previously-coalesced*, determine which lineages have coalesced (wrt a specific site). Update *previously-coalesced* accordingly. And then (not yet implemented) add the corresponding tree node to *tree-builder*."
;;   (let* ((all-descendants-of-lineage (find-descendants site edge))
;;          (coalescence-time (first edge))
;;          (lineage-pair
;;            (loop for x in *previously-coalesced*
;;                  with r = all-descendants-of-lineage ; r = remaining lineages
;;                  until (<= (length r) 1)
;;                  when (sort (intersection x r) #'<)
;;                    collect it into lineage-pair
;;                    and do (setf r (set-difference r (intersection r x)))
;;                  finally
;;                     (return (cond
;;                               ;; Case 1: both coalescing lineages are nodes
;;                               ((= (length lineage-pair) 2)
;;                                (progn
;;                                  (format t "~%Case 1.~%lineage-pair: ~a" lineage-pair)
;;                                  (setf *previously-coalesced*
;;                                        (cons
;;                                         (sort (copy-list
;;                                                (union (first lineage-pair)
;;                                                       (second lineage-pair)))
;;                                               #'<)
;;                                         *previously-coalesced*))
;;                                  lineage-pair))
;;                               ;; Case 2: one lineage is a node and the other is a leaf
;;                               ((and (= (length lineage-pair) 1) (= (length r) 1))
;;                                (progn
;;                                  (format t "~%Case 2.~%lineage-pair: ~a~%Leaf: ~a~%Node: ~a"
;;                                          lineage-pair r (first lineage-pair))
;;                                  (setf *previously-coalesced*
;;                                        (cons
;;                                         (sort (copy-list
;;                                                (union r (first lineage-pair)))
;;                                               #'<)
;;                                         *previously-coalesced*))
;;                                  (list (first lineage-pair) r)))
;;                               ;; Case 3: both lineages have already coalesced at
;;                               ;; a more recent time, so ignore this edge.
;;                               ((and (= (length lineage-pair) 1) (= (length r) 0)) nil)
;;                               ;; Case 4: both lineages are leaves
;;                               ((and (= (length lineage-pair) 0) (= (length r) 2))
;;                                (progn
;;                                  (format t "~%Case 4.~%lineage-pair: ~a" lineage-pair)
;;                                  (setf *previously-coalesced*
;;                                        (cons r *previously-coalesced*))
;;                                  r)))))))
;;     (format t "~%coalescing pair is: ~a" lineage-pair)
;;     (format t "~%all descendants are: ~a" all-descendants-of-lineage)
;;     ;;(format t "~%coalescence-time is: ~a" coalescence-time)
;;     (format t "~%updated *previously-coalesced*: ~a" *previously-coalesced*)))


;; ;; ;; Example
;; (defparameter test-node (make-node :population-name '(1 4) 
;;                                     :population-start-time 1.2
;;                                     :left-subtree (make-leaf :population-name 2
;;                                                              :population-end-time 1.2)
;;                                     :right-subtree (make-leaf :population-name 3
;;                                                               :population-end-time 1.2)))
;; (defparameter test-node2 (make-node :population-name '(2 3) 
;;                                     :population-start-time 2.2
;;                                     :left-subtree (make-leaf :population-name 2
;;                                                              :population-end-time 2.2)
;;                                     :right-subtree (make-leaf :population-name 3
;;                                                               :population-end-time 2.2)))

;; (defparameter *tree-builder* (acons '(1 4) (list test-node) nil))
;; (setf *tree-builder* (acons '(2 3) (list test-node2) *tree-builder*))

;; *tree-builder*
;; (((2 3)
;;   ((:POPULATION-NAME (2 3) :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE
;;     0 :PARENT-LABEL "parent" :POPULATION-START-TIME 2.2 :POPULATION-END-TIME
;;     NIL :NUMERIC-LABEL NIL)
;;    ((:POPULATION-NAME 2 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
;;      :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME 2.2
;;      :NUMERIC-LABEL NIL))
;;    ((:POPULATION-NAME 3 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
;;      :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME 2.2
;;      :NUMERIC-LABEL NIL))))
;;  ((1 4)
;;   ((:POPULATION-NAME (1 4) :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE
;;     0 :PARENT-LABEL "parent" :POPULATION-START-TIME 1.2 :POPULATION-END-TIME
;;     NIL :NUMERIC-LABEL NIL)
;;    ((:POPULATION-NAME 1 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
;;      :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME 1.2
;;      :NUMERIC-LABEL NIL))
;;    ((:POPULATION-NAME 4 :DIST-FROM-PARENT NIL :MUTATION-RATE 1 :RECOMB-RATE 0
;;      :PARENT-LABEL "parent" :POPULATION-START-TIME NIL :POPULATION-END-TIME 1.2
;;      :NUMERIC-LABEL NIL)))))

;; (assoc '(1 4) *tree-builder* :test #'equal)
;; (assoc '(2 3) *tree-builder* :test #'equal)
;; (nodep (second (assoc '(1 4) *tree-builder* :test #'equal)))



;; ;; Example
;; *list-of-breakpoints* = (402 335 766 754 973 186)

;; (mscr t* 1000)
;; (((4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000))))
;;  ((4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000)))
;;   (2.854556611726914d0 ((186 . 1000)) NIL NIL NIL)
;;   (3.003349457914102d0 ((335 . 1000)) NIL ((1 . 1000)) NIL)
;;   (3.0590606609388025d0 ((186 . 334)) ((754 . 1000)) NIL NIL)
;;   (3.2285537794507313d0 ((186 . 1000)) ((754 . 1000)) ((1 . 1000)) NIL)
;;   (3.362515819129609d0 ((186 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;;   (3.4149221114926416d0 ((1 . 185)) NIL NIL NIL)
;;   (3.715127784448965d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)))

;; (defparameter *sorted-edges*
;;   (process-mscr-output
;;    '(((4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000))))
;;      ((4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000)))
;;       (2.854556611726914d0 ((186 . 1000)) NIL NIL NIL)
;;       (3.003349457914102d0 ((335 . 1000)) NIL ((1 . 1000)) NIL)
;;       (3.0590606609388025d0 ((186 . 334)) ((754 . 1000)) NIL NIL)
;;       (3.2285537794507313d0 ((186 . 1000)) ((754 . 1000)) ((1 . 1000)) NIL)
;;       (3.362515819129609d0 ((186 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;;       (3.4149221114926416d0 ((1 . 185)) NIL NIL NIL)
;;       (3.715127784448965d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)))))

;; ((2.854556611726914d0 ((186 . 1000)) NIL NIL NIL)
;;  (3.003349457914102d0 ((335 . 1000)) NIL ((1 . 1000)) NIL)
;;  (3.0590606609388025d0 ((186 . 334)) ((754 . 1000)) NIL NIL)
;;  (3.2285537794507313d0 ((186 . 1000)) ((754 . 1000)) ((1 . 1000)) NIL)
;;  (3.362515819129609d0 ((186 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;;  (3.4149221114926416d0 ((1 . 185)) NIL NIL NIL)
;;  (3.715127784448965d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;;  (4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000))))

;; (progn (defparameter *previously-coalesced* nil)
;;        (loop for edge in *sorted-edges*
;;              do (progn
;;                   (format t "~%~%EDGE: ~a" edge)
;;                   (make-node-from-edge 401 edge))))


;; EDGE: (2.854556611726914d0 ((186 . 1000)) NIL NIL NIL)
;; coalescing pair is: NIL
;; all descendants are: (1)
;; updated *previously-coalesced*: NIL

;; EDGE: (3.003349457914102d0 ((335 . 1000)) NIL ((1 . 1000)) NIL)
;; Case 4.
;; lineage-pair: NIL
;; coalescing pair is: (1 3)
;; all descendants are: (1 3)
;; updated *previously-coalesced*: ((1 3))

;; EDGE: (3.0590606609388025d0 ((186 . 334)) ((754 . 1000)) NIL NIL)
;; coalescing pair is: NIL
;; all descendants are: NIL
;; updated *previously-coalesced*: ((1 3))

;; EDGE: (3.2285537794507313d0 ((186 . 1000)) ((754 . 1000)) ((1 . 1000)) NIL)
;; coalescing pair is: NIL
;; all descendants are: (1 3)
;; updated *previously-coalesced*: ((1 3))

;; EDGE: (3.362515819129609d0 ((186 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;; Case 2.
;; lineage-pair: ((1 3))
;; Leaf: (2)
;; Node: (1 3)
;; coalescing pair is: ((1 3) (2))
;; all descendants are: (1 2 3)
;; updated *previously-coalesced*: ((1 2 3) (1 3))

;; EDGE: (3.4149221114926416d0 ((1 . 185)) NIL NIL NIL)
;; coalescing pair is: NIL
;; all descendants are: NIL
;; updated *previously-coalesced*: ((1 2 3) (1 3))

;; EDGE: (3.715127784448965d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) NIL)
;; coalescing pair is: NIL
;; all descendants are: (1 2 3)
;; updated *previously-coalesced*: ((1 2 3) (1 3))

;; EDGE: (4.306453153051099d0 ((1 . 1000)) ((1 . 1000)) ((1 . 1000)) ((1 . 1000)))
;; Case 2.
;; lineage-pair: ((1 2 3))
;; Leaf: (4)
;; Node: (1 2 3)
;; coalescing pair is: ((1 2 3) (4))
;; all descendants are: (1 2 3 4)
;; updated *previously-coalesced*: ((1 2 3 4) (1 2 3) (1 3))









  (cond ((endp lineage-pair) ; two leaves are coalescing
         (setf *tree-builder*
               (acons all-descendants-of-lineage
                      (make-node :population-name (first all-descendants-of-lineage)
                                 :population-start-time coalescence-time
                                 :left-subtree
                                 (make-leaf :population-name (first all-descendants-of-lineage)
                                            :population-start-time "need to get this from species tree"
                                            :dist-from-parent "need to get this from species tree"
                                            :population-end-time coalescence-time)
                                 :right-subtree
                                 (make-leaf :population-name (second all-descendants-of-lineage)
                                            :population-start-time "need to get this from species tree"
                                            :dist-from-parent "need to get this from species tree"
                                            :population-end-time coalescence-time))
                      *previously-coalesced*)
               *previously-coalesced*
               (cons all-descendants-of-lineage *previously-coalesced*)))

        ((= 1 (length lineage-pair)) ; leaf is coalescing with a node
         (setf *tree-builder*
               (acons all-descendants-of-lineage
                      (make-node :population-name (first all-descendants-of-lineage)
                                 :population-start-time coalescence-time
                                 :left-subtree
                                 ; left off here
                                 (cdr (assoc something *tree-builder* :test #'equalp))
                                 (assoc lineage-pair
                                 :right-subtree
                                 (make-leaf :population-name (second all-descendants-of-lineage)
                                            :population-start-time "need to get this from species tree"
                                            :dist-from-parent "need to get this from species tree"
                                            :population-end-time coalescence-time))
                      *previously-coalesced*)
               *previously-coalesced*
               (cons all-descendants-of-lineage *previously-coalesced*)))
        
         ((= 2 (length lineage-pair)) "case 2") ; two nodes are coalescing
         (t "if you got here something went wrong"))
        )

