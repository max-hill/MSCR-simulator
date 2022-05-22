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

;; Part 1. Constructors
;;
;; There are two kinds of binary trees: (1) nodes and (2) leaves. 
;;
;; 1. A *leaf* is a list with two components: a string (i.e. the label) followed
;;    by a list of three elements (1: a number representing the distance from
;;    the parent, 2: a number, representing the mutation rate on that edge, and
;;    3: a number, representing the recombination rate on that edge). For
;;    example, ("a" (2.3 .1 1.2))
;;
;; 2. A *node* is a list with three components: a leaf, a left subtree and a
;;    right subtree: (leaf left-subtree right-subtree). For example:
;;    (("ab" 1) ("a" 2.3) ("b" 2.7))


(defun make-leaf (leaf-label distance-from-parent mutation-rate recombination-rate)
  "Create a leaf. Leaf-label should be a string. Distance-from-parent should be
a number."
  (list leaf-label
        (list distance-from-parent
        mutation-rate
        recombination-rate)))

(defun make-node (node-label distance-from-parent mutation-rate
recombination-rate left-subtree right-subtree)
  "Create a node. Node-label should be a string, distance-from-parent should be
a number, and the subtrees should be trees (either nodes or leafs, but not nil."
  (list (make-leaf node-label
                   distance-from-parent
                   mutation-rate
                   recombination-rate)
        left-subtree
        right-subtree))

;; Example 1. A rooted three-leaf tree: (deprecated)
;; (make-node "root" 0
;;            (make-node "AB" 1
;;                       (make-leaf "A" 1)
;;                       (make-leaf "B" 1))
;;            (make-leaf "C" 2))
;;
;; Example 2. A balanced quartet: (deprecated)
;; (make-node "root" 0
;;            (make-node "AB" 1
;;                       (make-leaf "A" .1)
;;                       (make-leaf "B" .1))
;;            (make-node "CD" 1
;;                       (make-leaf "C" .1)
;;                       (make-leaf "D" .1)))
;;
;; Example 3. An unbalanced quartet:
;; (make-node "root" 0 3 0 
;;            (make-node "ABC" 1 2 3 
;;                       (make-node "AB" .1 2 3
;;                                  (make-leaf "A" .1 1 4)
;;                                  (make-leaf "B" .1 1 1))
;;                       (make-leaf "C" .1 1.3 0))
;;            (make-leaf "D" 4 1.1 0))



;; Part 2. Recognizers

(defun leafp (tree)
  "Return t if tree is a leaf and nil otherwise."
  (if (and (listp tree)
           (= 2 (length tree))
           (stringp (first tree))
           (listp (second tree))
           (= 3 (length (second tree))))
      t
      nil))

(defun nodep (tree)
  "Return t if tree is a node and nil otherwise."
  (if (and (listp tree)
           (= 3 (length tree)))
      t
      nil))

;; Part 3. Selectors

(defun leaf-name (leaf)
  "Return name of leaf."
  (first leaf))


(defun leaf-dist-from-parent (leaf)
  "Return distance of leaf from its parent."
  (first (second leaf)))

(defun node-name (node)
  "Return the name (i.e. label) of node."
  (first (first node)))

(defun get-leaf-number (tree))
  
(defun fast-bin-tree-preorder (tree)
  "A tail-recursive version of bin-tree-preorder. Source:
https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html"
  (preorder-aux tree nil))

(defun list-leaves-aux (tree x)
  "recursively build a list of leaf names of tree. the accumulator variable is x. Source: https://www2.cs.sfu.ca/CourseCentral/310/pwfong/Lisp/3/tutorial3.html"
  (if (leafp tree)
      (cons (leaf-name tree) x)
      (cons (list-leaves-aux (right-subtree tree)
                             (list-leaves-aux (left-subtree tree) x)))))




  
(defun node-dist-from-parent (node)
  "Return the distance of node from its parent."
  (first (second (first node))))

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
  (if (nodep tree)
      (node-name tree)
      (leaf-name tree)))

(defun get-vertex-dist-from-parent (tree)
  (if (nodep tree)
      (node-dist-from-parent tree)
      (leaf-dist-from-parent tree)))


;; Part 4. Simulator: This part simulates the JC process on a (general) tree.
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

;; We'll need a way to perform nucleotide substitutions when they are determined
;; to occur. This functionality is provided by the next function.

(defun substitute-nucleotide (current-nucleotide)
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
         (edge-length (get-vertex-dist-from-parent tree))
         (substitution-probability
           (* .75 (- 1 (exp (* (/ -4 3) *mutation-rate* edge-length)))))
         (new-nucleotide (if (< (random 1d0) substitution-probability)
                             (substitute-nucleotide parent-nucleotide)
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
  "Tail-recursive implmentation of JC69 on the tree. Starts with random
nucleotide at the root."
  (let ((root-state (draw-random-nucleotide)))
    (progn (format t "root state: ~a ~%" root-state)
           (evolve-down-tree-aux tree root-state nil))))

(defun evolve-down-tree-aux (tree parent-nucleotide x)
  "x is the accumulator variable"
  (if (leafp tree)
      (cons (evolve-down-edge tree parent-nucleotide) x)
      (let ((new-parent-nucleotide
              (second (evolve-down-edge tree parent-nucleotide)))
            (left (left-subtree tree))
            (right (right-subtree tree)))
        (evolve-down-tree-aux
         left
         new-parent-nucleotide
         (evolve-down-tree-aux right new-parent-nucleotide x)))))

;; It works! Example usage:
;; (evolve-down-tree (make-node "root" 0
;;                              (make-node "AB" 1
;;                                         (make-leaf "A" .1)
;;                                         (make-leaf "B" .1))
;;                              (make-leaf "C" 1)))



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


(((0 (1 2 3 4 5) NIL NIL)))


Make sample
(make-list 5)
(list 0 (interval 1 5))


;; the following function allows us to avoid the strange backticks used in
;; simulate-three-species
(defun make-leaf-sample (i n k)
  "makes the i-th initial sample for an n-taxa gene tree when sequences are
length k"
  (list (list (cons 0 (append (make-list (- i 1))
                              (cons (interval 1 k)
                                    (make-list (- n i))))))))


(defun simulate-three-species (τ_ab τ_abc τ_max ρ_a ρ_b ρ_c ρ_ab ρ_abc number-of-base-pairs)
  "Construct an ARG on a 3-taxa species tree with given parameters"
  (let* ((sample-A (make-leaf-sample 1 3 number-of-base-pairs))
	 (sample-B (make-leaf-sample 2 3 number-of-base-pairs))
	 (sample-C (make-leaf-sample 3 3 number-of-base-pairs))
	 (edges-from-A (arg-builder ρ_a 0 τ_ab number-of-base-pairs sample-A))
	 (edges-from-B (arg-builder ρ_b 0 τ_ab number-of-base-pairs sample-B))
	 (edges-from-C (arg-builder ρ_c 0 τ_abc number-of-base-pairs sample-C))
	 (edges-from-AB (arg-builder ρ_ab τ_ab τ_abc number-of-base-pairs
				     (list (union (first edges-from-A)
						  (first edges-from-B))
					   (union (second edges-from-A)
						  (second edges-from-B))))))
    (arg-builder ρ_abc τ_abc τ_max number-of-base-pairs
		 (list (union (first edges-from-AB)
			      (first edges-from-C))
		       (union (second edges-from-AB)
			      (second edges-from-C)))
		 t)))


(defun merge-output-edge-sets (s1 s2)
  "combine active lineages from two daughter populations for entry into their
parent population. Note that if one of the edge-sets is nil, this will just
output the other edge set."
  (list (union (first s1) (first s2))))

(defun mscr-general-simulator (i n number-of-base-pairs tree edge-sets)
  "n: number of leaves of initial tree. number-of-base-pairs: sequence length, in nucleotides. The
input i tracks the leaf number: its initial value is 0 and it increments i by 1
each time a right subtree is recursively evaluated. initial value of edge-sets
should be nil, i think."
  (if (leafp tree)
      (make-leaf-sample leaf-number-tracker n number-of-base-pairs)
      (arg-builder n number-of-base-pairs i tree
                   (merge-output-edge-sets (mscr-general-simulator n number-of-base-pairs i
                                                                   (left-subtree tree) edge-sets)
                                           (mscr-general-simulator n number-of-base-pairs (+ i 1)
                                                                   (right-subtree tree) edge-sets)))))
      
(let* ((k (number-of-active-lineages edge-sets))
       (coales-rate (* .5 k (- k 1)))
       (recomb-rate (* .5 k ρ))
       (total-rate (+ coales-rate recomb-rate))
       (t_1 (+ t_0 (draw-exponential total-rate))))
  (if (or (> t_1 t_end) (and stop-at-mrca (= k 1)))
      edge-sets
      (arg-builder ρ t_1 t_end number-of-base-pairs
		   (if (< (random 1d0) (/ recomb-rate total-rate))
		       (implement-recombination t_1 edge-sets number-of-base-pairs)
		       (implement-coalescence t_1 edge-sets))
		   stop-at-mrca)))

(defun number-of-active-lineages (edge-sets)
  "Determine how many lineages are 'active' in a population (a lineage is active
if it is a candidate for coalescence. Recall that when a lineage undergoes
recombination, it is rendered inactive, but two new active lineages are created.
In general, active lineages are stored in P and inactive lineages are stored in
Q. I can't remember off the top of my head but there might be an exception to
this rule.)"
  (length (first edge-sets)))
