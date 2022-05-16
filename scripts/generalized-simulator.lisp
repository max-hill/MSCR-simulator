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
;;    by a number (i.e. the distance from the parent). For example, ("a" 2.3)
;;
;; 2. A *node* is a list with three components: a leaf, a left subtree and a
;;    right subtree: (leaf left-subtree right-subtree). For example:
;;    (("ab" 1) ("a" 2.3) ("b" 2.7))


(defun make-leaf (leaf-label distance-from-parent)
  "Create a leaf. Leaf-label should be a string. Distance-from-parent should be
a number."
  (list leaf-label distance-from-parent))

(defun make-node (node-label distance-from-parent left-subtree right-subtree)
  "Create a node. Node-label should be a string, distance-from-parent should be
a number, and the subtrees should be trees (either nodes or leafs, but not nil."
  (list (make-leaf node-label distance-from-parent) left-subtree right-subtree))

;; For example,
;; (make-node "root" 0
;;            (make-node "AB" 1
;;                       (make-leaf "A" 1)
;;                       (make-leaf "B" 1))
;;            (make-leaf "C" 2))


;; Part 2. Recognizers

(defun leafp (tree)
  "Return t if tree is a leaf and nil otherwise."
  (if (and (listp tree)
           (= 2 (length tree))
           (stringp (first tree))
           (numberp (second tree)))
      t
      nil))

(defun nodep (tree)
  "Return t if tree is a node and nil otherwise."
  (if (and (listp tree)
           (= 3 (length tree)))
      t
      nil))

(defun leaf-name (leaf)
  "Return name of leaf."
  (first leaf))


;; Part 3. Selectors

(defun leaf-dist-from-parent (leaf)
  "Return distance of leaf from its parent."
  (second leaf))

(defun node-name (node)
  "Return the name (i.e. label) of node."
  (first (first node)))


(defun node-dist-from-parent (node)
  "Return the distance of node from its parent."
  (second (first node)))

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
;; 1. Implement a newick -> lisp tree converter
;; 2. Actually implement recombination...
;; 3. Write an msr. 
