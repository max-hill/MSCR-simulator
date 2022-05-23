;;;; generalized-simulator.lisp --- My attempt to implement MSCR on general
;;;;                                trees. Currently under construction.
;; 
;; Author: max hill 
;; (Last updated 2022-05-15)

;; DESCRIPTION: Here our goal is to implement a simulator which takes as input a
;; binary tree T (ideally in newick tree format, with branch lengths) and output
;; an MSA in the form of an alignment file generated according to the MSCR
;; process on T.
;;
;; We'll start by following the ideas presented at
;; https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html
;;
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

(defun make-leaf (leaf-label distance-from-parent mutation-rate recombination-rate)
  "Create a leaf. Leaf-label should be a string. Distance-from-parent should be
a number."
  `(,leaf-label
    (,distance-from-parent
     ,mutation-rate
     ,recombination-rate)))

(defun make-node (node-label distance-from-parent mutation-rate
                  recombination-rate left-subtree right-subtree)
  "Create a node. Node-label should be a string, distance-from-parent should be
a number, and the subtrees should be trees (either nodes or leafs, but not nil."
  `((,node-label (,distance-from-parent ,mutation-rate ,recombination-rate))
    ,left-subtree
    ,right-subtree))

;;______________________________________________________________________________
;;
;; Part 2. Recognizers
;;______________________________________________________________________________

(defun leafp (tree)
  "Return t if tree is a leaf and nil otherwise. Does not distinguish between
leaves with/without vertex times."
  (if (and (listp tree)
           (= 2 (length tree))
           (stringp (first tree))
           (listp (second tree))
           (or (= 4 (length (second tree)))
               (= 3 (length (second tree)))))
      t
      nil))

(defun nodep (tree)
  "Return t if tree is a node and nil otherwise. Does not distinguish between
  leaves with/without vertex times."
  (if (and (listp tree)
           (= 3 (length tree))
           (or (nodep (second tree))
               (leafp (second tree)))
           (or (nodep (third tree))
               (leafp (third tree))))
      t
      nil))

;;______________________________________________________________________________
;;
;; Part 3. Add vertex times to the trees
;;______________________________________________________________________________

(defun get-edge-parameters (tree)
  "Return the list of parameters corresponding to the edge above the vertex"
  (if (leafp tree)
      (second tree)
      (second (first tree))))

(defun add-vertex-time (tree parent-age)
  "Adds the vertex time to a leaf object."
  (list (get-vertex-name tree)
        (cons (- parent-age (get-dist-from-parent tree))
              (get-edge-parameters tree))))

(defun add-vertex-times-to-tree-aux (tree parent-age)
  "Recursive auxilliary function for the function add-vertex-times-to-tree."
  (if (leafp tree)
      (add-vertex-time tree parent-age)
      (list (add-vertex-time tree parent-age)
            (add-vertex-times-to-tree-aux (left-subtree tree)
                                          (- parent-age
                                             (get-dist-from-parent tree)))
            (add-vertex-times-to-tree-aux (right-subtree tree)
                                          (- parent-age
                                             (get-dist-from-parent tree))))))

(defun get-total-tree-height (tree)
  "Returns the total height of the tree, i.e. from the top of the root edge to
the maximally distant leaf."
  (get-tree-height-aux tree 0))

(defun get-tree-height-aux (tree accumulator)
  "Returns the total height of the tree, i.e. the length of the maximumal path
from the vertex to one of its leaves"
  (if (leafp tree)
      (+ accumulator (get-dist-from-parent tree))
      (max (get-tree-height-aux (right-subtree tree)
                                (+ accumulator (get-dist-from-parent tree)))
           (get-tree-height-aux (left-subtree tree)
                                (+ accumulator (get-dist-from-parent tree))))))

(defun add-vertex-times-to-tree (tree)
  "Input: a tree constructed with the functions 'make-node' and 'make-leaf'.
Output: the same tree, but with labelled vertex times."
  (add-vertex-times-to-tree-aux tree (get-total-tree-height tree)))

        
;;______________________________________________________________________________
;;
;; Part 4. Tree examples
;;______________________________________________________________________________
;; From here on we will assume that all trees have vertex times, i.e. were
;; processed with the function 'add-vertex-times-to-tree'

;; Example 1. A leaf
(defparameter *l1*
  (make-leaf "A" .1 1 4))
(defparameter *lv1*
  (add-vertex-times-to-tree *l1*))

;; Example 2. A node (unbalanced quartet):
(defparameter *n1*
  (make-node "root" 0 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv1* (add-vertex-times-to-tree *n1*))

;; Example 3. A node (unbalanced quartet):
(defparameter *n2*
  (make-node "root" 999 3 0
             (make-node "ABC" 1 2 3 
                        (make-node "AB" .1 2 3
                                   (make-leaf "A" .1 1 4)
                                   (make-leaf "B" .1 1 1))
                        (make-leaf "C" .1 1.3 0))
             (make-leaf "D" 4 1.1 0)))
(defparameter *nv2* (add-vertex-times-to-tree *n2*))

;; Example 4. A node (balanced quartet):
(defparameter *n3*
  (make-node "root" 999 .1 .2
             (make-node "ab" 2 .1 .2
                        (make-leaf "a" 2 .1 .2)
                        (make-leaf "b" 2 .1 .2))
             (make-node "cd" 1 .1 .2
                        (make-leaf "c" 3 .1 .2)
                        (make-leaf "d" 2 .1 .2))))
(defparameter *nv3* (add-vertex-times-to-tree *n3*))

;; Example 5. Three taxa example
(defparameter *3-taxa-example* '(("abc" (4 999 1 0))
                                 (("ab" (2 2 0.1 0.2))
                                  ("a" (0 2 0.1 0.2))
                                  ("b" (0 2 0.1 0.2)))
                                 ("c" (0 4 1 0))))

;;______________________________________________________________________________
;;
;; Part 5. More Selectors
;;______________________________________________________________________________
;; all these assume that the tree contains vertex times

;; Leaf-specific selectors
(defun leaf-name (leaf)
  "Return name of leaf."
  (first leaf))
(defun leaf-dist-from-parent (leaf)
  "Return distance of leaf from its parent."
  (second (second leaf)))
(defun leaf-mutation-rate (leaf)
  "Return the mutation rate on the edge from leaf to its parent."
  (third (second leaf)))
(defun leaf-recombination-rate (leaf)
  "Return the recombination rate on the edge from leaf to its parent."
  (fourth(second leaf)))

;; Node-specific selectors
(defun node-name (node)
  "Return the name (i.e. label) of node."
  (first (first node)))
(defun node-dist-from-parent (node)
  "Return the length of the edge connecting node to its parent."
  (leaf-dist-from-parent (first node)))
(defun node-mutation-rate (node)
  "Return the mutation rate of the edge connecting node to its parent."
  (leaf-mutation-rate (first node)))
(defun node-recombination-rate (node)
  "Return the recombination rate of the edge connecting node to its parent."
  (leaf-recombination-rate (first node)))

;; General selectors
(defun left-subtree (tree)
  "Return the left subtree of tree. If tree is a leaf, return nil."
  (if (nodep tree)
      (second tree)
      nil))
(defun right-subtree (tree)
  "Return the right subtree of tree. If tree is a leaf, return nil."
  (if (nodep tree)
      (third tree)
      nil))
(defun get-vertex-name (tree)
  "Return the label of vertex. Works for both leaves and internal nodes."
  (if (nodep tree)
      (node-name tree)
      (leaf-name tree)))
(defun get-dist-from-parent (tree)
  "Return the length of the edge connecting vertex to its parent. Works for both
leaves and internal nodes."
  (if (nodep tree)
      (node-dist-from-parent tree)
      (leaf-dist-from-parent tree)))
(defun get-mutation-rate (tree)
  "Return the mutation rate along the edge connecting vertex to its parent.
Works for both leaves and internal nodes."
  (if (nodep tree)
      (node-mutation-rate tree)
      (leaf-mutation-rate tree)))
(defun get-recomb-rate (tree)
  "Return the recombination rate along the edge connecting vertex to its parent.
Works for both leaves and internal nodes."
  (if (nodep tree)
      (node-recombination-rate tree)
      (leaf-recombination-rate tree)))
(defun get-vertex-time (tree)
  "Return the vertex time, i.e. the age of the vertex."
  (if (leafp tree)
      (first (second tree))
      (first (second (first tree)))))

;;______________________________________________________________________________
;;
;; Part 6. Mutation Simulator.
;;______________________________________________________________________________
;; This part simulates the JC process on a (general) tree.
;;
;; For simplicity, we will implement nucleotides here as the numbers 0,1,2,3
;; rather than A,T,C,G. Everything will assume JC69 process.
;;
;; We'll just us a global mutation rate for now. This is the rate at which a
;; nucleotide leaves is current state according to the JC69 process. It is
;; constant for the whole species tree.

(defparameter *mutation-rate* 1)

;; The next function is used for generating a nucleotide at the root of the tree.

(defun draw-random-nucleotide ()
  "Return a random number from the set {0,1,2,3}."
  (random 4))

;; The next (highly inefficient) function creates an interval of integers.

(defun interval (a b)
  "Construct an 'interval' of integers from a to b, inclusive. Note that a and b
must both be integers."
  (loop for i from a to b collecting i))
 
;; We'll need a way to perform nucleotide substitutions when they are determined
;; to occur. This functionality is provided by the next function.

(defun implement-substitution (current-nucleotide)
  "Implements a substitution to a *different* nucleotide. Specifically, return a
random number from the set {0,1,2,3}\{current-nucleotide}. It does this by
adding 1d3 to the current nucleotide and then reducing the result modulo 4."
  (declare (integer current-nucleotide))
  (mod (+ current-nucleotide 1 (random 3)) 4))

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
         (edge-length (get-dist-from-parent tree))
         (mutation-rate (get-mutation-rate tree))
         (substitution-probability
           (* .75 (- 1 (exp (* (/ -4 3) mutation-rate edge-length)))))
         (new-nucleotide (if (< (random 1d0) substitution-probability)
                             (implement-substitution parent-nucleotide)
                             parent-nucleotide)))
    (progn
      (if (nodep tree)
          (unless (equal vertex-name "root")
            (format t "ancestor ~a: ~a ~%" vertex-name new-nucleotide))
          (format t "leaf ~a: ~a ~%" vertex-name new-nucleotide))
      (list vertex-name new-nucleotide))))

;; The next two functions are the workhorses of our simulator. They use the
;; above functions to execute a tail-recursive implmentation of the JC69 process
;; on the input tree.

(defun evolve-down-tree (tree)
  "Recursive implmentation of JC69 on the tree. Starts with random nucleotide at
the tip of the edge above the root. To make this start exactly at the root (by
which I mean the mrca of the samples) you should set the 'dist-from-parent' in
the root vertex to zero."
  (let ((root-state (draw-random-nucleotide)))
    (progn (format t "root state: ~a ~%" root-state)
           (evolve-down-tree-aux tree root-state nil))))

(defun evolve-down-tree-aux (tree parent-nucleotide x)
  "x is the accumulator variable use to construct the output list"
  (if (leafp tree)
      (cons (evolve-down-edge tree parent-nucleotide) x)
      (let ((new-parent-nucleotide (second (evolve-down-edge tree parent-nucleotide))))
        (evolve-down-tree-aux (left-subtree tree)
                              new-parent-nucleotide
                              (evolve-down-tree-aux (right-subtree tree)
                                                    new-parent-nucleotide
                                                    x)))))

;; It works! Example usage:
;; (evolve-down-tree *example-node*)


;; Remains to do:
;;
;; 1. Implement a newick -> lisp tree converter
;;
;; 2. Implement recursion. To do this, we need to do two tasks. First, we need
;;    to adapt the main arg-simulator function arg-builder to a similar
;;    recursive structure as we have done here in evolve-down-tree. This new
;;    recursive arg-builder function will output a collection called
;;    'output-edges' of coalescent times of all the pairs. Second, to generate a
;;    marginal gene tree for each of the k sites, we can run
;;    (compute-marginal-tmrcas i j output-edges k) for all distinct pairs of
;;    leaves i and j. Then we'll have a big matrix of coalescent times and we
;;    can just lookup the times to create the marginal gene trees. Or better, we
;;    can write a new function that takes as input the set output-edges, a site
;;    i∈{1,...k}, and outputs a gene tree in the appropriate form for
;;    evolve-down-tree. Then run evolve down tree and save the result as a
;;    column in the MSA.
;;
;; 3. Implement the inference methods on the new MSAs. In particular, we want to
;;    implement the four-point method. We could also try estimating branch
;;    lengths, as this would be a key test of one of our paper's conclusions.


;; the following function allows us to avoid the strange backticks used in
;; simulate-three-species
(defun make-leaf-sample (i n-total-leaf-number k-sequence-length)
  "makes the i-th initial sample for a phylogenetic tree with a number of taxa
equal to n-total-leaf-number, and when sequences are length k-sequence-length"
  (list (list (cons 0 (append (make-list (- i 1))
                              (cons (interval 1 k-sequence-length)
                                    (make-list (- n-total-leaf-number i))))))))


(defun simulate-three-species (τ_ab τ_abc τ_max ρ_a ρ_b ρ_c ρ_ab ρ_abc k-sequence-length)
  "Construct an ARG on a 3-taxa species tree with given parameters"
  (let* ((sample-A (make-leaf-sample 1 3 k-sequence-length))
	 (sample-B (make-leaf-sample 2 3 k-sequence-length))
	 (sample-C (make-leaf-sample 3 3 k-sequence-length))
	 (edges-from-A (arg-builder ρ_a 0 τ_ab k-sequence-length sample-A))
	 (edges-from-B (arg-builder ρ_b 0 τ_ab k-sequence-length sample-B))
	 (edges-from-C (arg-builder ρ_c 0 τ_abc k-sequence-length sample-C))
	 (edges-from-AB (arg-builder ρ_ab τ_ab τ_abc k-sequence-length
				     (list (union (first edges-from-A)
						  (first edges-from-B))
					   (union (second edges-from-A)
						  (second edges-from-B))))))
    (arg-builder ρ_abc τ_abc τ_max k-sequence-length
		 (list (union (first edges-from-AB)
			      (first edges-from-C))
		       (union (second edges-from-AB)
			      (second edges-from-C)))
		 t)))


(defun merge-output-edge-sets (child-1 child-2)
  "Combine the output edge sets child-1=(P_1,Q_1) and child-2=(P_2,Q_2) from two
daughter populations for entry into their parent population. Output a pair
of (P,Q) of active and inactive edges in the appropriate format for use in
arg-builder. Note that if one of the edge-sets is nil, this will just output the
other edge set."
  (list (union (first child-1) (first child-2))
        (union (second child-1) (second child-2))))


(defun get-population-start-time (tree)
  (get-vertex-time tree))
(defun get-population-end-time (tree)
  (+ (get-vertex-time tree)
     (get-dist-from-parent tree)))

;; this is the main simulator. It is currently running without error, but the
;; output is not correct. It looks like some of the subroutines borrowed from the earlier simulator need to be modified.
(defun mscr (n-total-leaf-number leaf-number-tracker k-sequence-length tree edge-sets)
  "Input comments: the variable 'n-total-leaf-number' is number of leaves on the
initial (i.e. full) tree, 'k-sequence-length' is just the length of the
sequences in base pairs. The input 'leaf-number-tracker' tracks the leaf number:
its initial value is 1? and it increments leaf-number-tracker by 1 each time a
right subtree is recursively evaluated. The initial value of edge-sets should be
nil, I think."
  (build-single-population-arg (get-recomb-rate tree)           ; ρ
                               (get-population-start-time tree) ; t_0
                               (get-population-end-time tree)   ; t_end
                               k-sequence-length
                               (if (leafp tree)
                                   (make-leaf-sample leaf-number-tracker
                                                     n-total-leaf-number
                                                     k-sequence-length)
                                   (merge-output-edge-sets
                                    (mscr n-total-leaf-number
                                          leaf-number-tracker
                                          k-sequence-length
                                          (left-subtree tree)
                                          edge-sets) ; is this right?
                                    (mscr n-total-leaf-number
                                          (+ leaf-number-tracker 1)
                                          k-sequence-length
                                          (right-subtree tree)
                                          edge-sets)))))

; this is the auxillary function to the main simulator. It is analogous to
; 'arg-builder' in the previous work
(defun build-single-population-arg (ρ t_0 t_end k-sequence-length edge-sets &optional (stop-at-mrca nil))
  "Build an ancestral recombination graph in a single population with fixed
start time t_0 and end time t_end."
  (let* ((number-of-active-lineages (get-number-of-active-lineages edge-sets))
	 (coales-rate (* .5 number-of-active-lineages (- number-of-active-lineages 1)))
	 (recomb-rate (* .5 number-of-active-lineages ρ))
	 (total-rate (+ coales-rate recomb-rate))
	 (t_1 (+ t_0 (draw-exponential total-rate))))
    (if (or (> t_1 t_end) (and stop-at-mrca (= number-of-active-lineages 1)))
	edge-sets
	(build-single-population-arg
         ρ
         t_1
         t_end
         k-sequence-length
	 (if (< (random 1d0) (/ recomb-rate total-rate))
	     (implement-recombination t_1 edge-sets k-sequence-length)
	     (implement-coalescence t_1 edge-sets))
	 stop-at-mrca))))



(defun get-number-of-active-lineages (edge-sets)
  "Determine how many lineages are 'active' in a population (a lineage is active
if it is a candidate for coalescence. Recall that when a lineage undergoes
recombination, it is rendered inactive, but two new active lineages are created.
In general, active lineages are stored in P and inactive lineages are stored in
Q. I can't remember off the top of my head but there might be an exception to
this rule.)"
  (length (first edge-sets)))



;; (defun fast-bin-tree-preorder (tree)
;;   "A tail-recursive version of bin-tree-preorder. Source:
;; https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html"
;;   (preorder-aux tree nil))

;; (defun list-leaves-aux (tree x)
;;   "recursively build a list of leaf names of tree. the accumulator variable is x. Source: https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html"
;;   (if (leafp tree)
;;       (cons (leaf-name tree) x)
;;       (cons (list-leaves-aux (right-subtree tree)
;;                              (list-leaves-aux (left-subtree tree) x)))))
