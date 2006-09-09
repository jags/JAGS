#ifndef INDEX_H_
#define INDEX_H_

/**
 * An Index behaves like a one-dimensional array of long integers, with
 * boundary checking.  Its principal uses are indexing the elements of an
 * SArray, and representing its dimension.
 *
 * @short array of long integers
 */
class Index {
  unsigned int _size;
  long *_index;
 public:
  /**
   * Default constructor. Creates a NULL Index of size 0
   */
  Index();
  /**
   * Constructs an index of given size. All elements of the index
   * are set to 1.
   */
  Index(unsigned int size);
  /**
   * Copy constructor
   */
  Index(Index const &index);
  virtual ~Index();
  /**
   * Assignment operator
   */
  Index &operator=(Index const &rhs);
  /**
   * Number of elements in the index
   */
  unsigned int size() const;
  /**
   * Subscripting operator. Allows an index to behave like an array.
   * @param i index to access (starting from 0). 
   * @exception out_of_range
   */
  long &operator[](unsigned int i);
  long operator[](unsigned int i) const;
  /**
   * Test for equality of two indices.
   */
  bool operator==(Index const &other) const;
  /**
   * Inequality operator.
   */
  bool operator!=(Index const &other) const;
  /**
   * Lexicographic ordering
   */
  bool operator<(Index const &other) const;
};

#endif /* INDEX_H_ */
