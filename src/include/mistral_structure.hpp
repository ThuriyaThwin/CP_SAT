
/*
  Mistral 2.0 is a constraint satisfaction and optimisation library
  Copyright (C) 2009  Emmanuel Hebrard
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  The author can be contacted electronically at emmanuel.hebrard@gmail.com.
*/


/*! \file mistral_structure.hpp
    \brief Header file for the basic data-structures.
*/


#include <iostream>
#include <iomanip>
#include <stdlib.h> 
#include <sstream> 


#ifndef __STRUCTURE_HPP
#define __STRUCTURE_HPP


#include <mistral_global.hpp>


namespace Mistral {


  /**********************************************
   * Vector
   **********************************************/
  /*! \class Vector
    \brief Simple vector class     
  */
  
  template <class DATA_TYPE>
  class Vector {
  public:

    /*!@name Parameters*/
    //@{
    DATA_TYPE* stack_;
    unsigned int capacity;
    unsigned int size;
    //@}

    /*!@name Constructor*/
    //@{
    Vector()
    {
      capacity = 0;
      size = 0;
      stack_ = NULL;
    }

    Vector(const Vector<DATA_TYPE>& s)
    {
      capacity = s.size;
      stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      for(size = 0; size<capacity; ++size)
	stack_[size] = s[size];
    }

    Vector(const int n)
    {
      capacity = n;
      size = 0;
      if( capacity )
	stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
      else stack_ = NULL;
    }
    //@}

    /*!@name Destructor*/
    //@{  
    virtual ~Vector()
    {
      free( stack_ );
    }
    //@}

    /*!@name Initialisation*/
    //@{
    void initialise(const unsigned int s, const unsigned int c)
    {
      size = s;
      capacity = c;
      stack_ = (DATA_TYPE*) malloc(capacity*sizeof(DATA_TYPE));
    }

    void extendStack()
    {
      capacity = ((capacity+1) << 1);
      DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
      stack_ = new_stack;
    }

    void extendStack( const unsigned int l )
    {
      capacity += l;
      DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
      stack_ = new_stack;
    }

    void resize( const unsigned int l )
    {
      if( capacity < l ) {
	capacity = l;
	DATA_TYPE* new_stack = (DATA_TYPE*) realloc(stack_, capacity*sizeof(DATA_TYPE));
	stack_ = new_stack;
      }
      size = l;
    }
    //@}

    /*!@name Accessors*/
    //@{
    inline int empty() const
    {
      return !size;
    }
  
    inline void add(DATA_TYPE x)
    {
      if( capacity == size ) 
	extendStack();
      stack_[size++] = x;
    }

    inline void secure(const DATA_TYPE x) 
    {
      if( capacity == size ) 
	extendStack();
    }

    inline void fast_add(DATA_TYPE x)
    {
      stack_[size++] = x;
    }

    inline void push_back(DATA_TYPE x)
    {
      if( capacity == size )
	extendStack();
      stack_[size++] = x;
    }

    inline DATA_TYPE popUntil(const unsigned int level)
    {
      size = level;
      return stack_[size];
    }

    inline DATA_TYPE pop()
    {
      return stack_[--size];
    }

    inline void pop(DATA_TYPE& x)
    {
      x = stack_[--size];
    }

    inline void clear()
    {
      size = 0;
    }

//     inline void remove(DATA_TYPE& elt)
//     {
//       unsigned int j=size;
//       while(j && stack_[--j] != elt);
//       stack_[j] = stack_[--size];
//     }
  
//     inline void erase(const unsigned int i)
//     {  
//       stack_[i] = stack_[--size];
//     }

    inline void remove(const unsigned int i)
    {  
      stack_[i] = stack_[--size];
    }

    inline void setBack(const DATA_TYPE& x, const int k=1)
    {
      stack_[size-k] = x;
    }

    inline DATA_TYPE& back(const int k=1)
    {
      return stack_[size-k];
    }

    inline const DATA_TYPE back(const int k=1) const
    {
      return stack_[size-k];
    }

    inline DATA_TYPE& operator[](const int i)
    {
      return stack_[i];
    }

    inline const DATA_TYPE operator[](const int i) const
    {
      return stack_[i];
    }

    inline Vector< DATA_TYPE >& operator=(const Vector< DATA_TYPE >& x)
    {
      initialise(0, x.capacity);
      for(unsigned int i=0; i<x.size; ++i)
	add(x[i]);
      return *this;
    }
    //@}

    /*!@name Printing*/
    //@{
    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << stack_[0] ;
      for(unsigned int i=1; i<size; ++i)
	os << " " << stack_[i];
      os << "]";
      return os;
    }
    //@}

  };


  /**********************************************
   * Array
   **********************************************/    
  /*! \class Array
    \brief Simple array class 
  */
  
  template <class DATA_TYPE>
  class Array {

  public: 

    /*!@name Parameters*/
    //@{
    int size;
    DATA_TYPE data[0];
    //@}

    Array(const Vector<DATA_TYPE>& ps) 
    {
      size = ps.size;
      for (int i=0; i<ps.size; ++i) 
	data[i] = ps[i];
    }

    static Array<DATA_TYPE>* Array_new(const Vector<DATA_TYPE>& ps)
    {
      void* mem = malloc(sizeof(Array<DATA_TYPE>) + sizeof(DATA_TYPE)*(ps.size));
      return new (mem) Array<DATA_TYPE>(ps); 
    }

    inline DATA_TYPE& operator [] (int i) { return data[i]; }
    inline DATA_TYPE operator [] (int i) const { return data[i]; }

    std::ostream& display(std::ostream& os) const {
      os << "[";
      if(size) os << data[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << data[i];
      os << "]";
      return os;
    }

  };


  /**********************************************
   * MultiSet
   **********************************************/
  /// Sparse multi-set representation

  class MultiSet 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    int *values_; // contiguous list of values (repeated)
    unsigned int **index_; // index[v][o] is the index of occurrence 'o' of value 'v' in 'values_'
    unsigned int *occurrence_; // self explanatory
    
    /// current max capacity
    unsigned int capacity;
    /// current size (total number of occurrences)
    unsigned int size;

    /// minimum and maximum value 
    int lb_;
    int ub_;
    //@}

    /*!@name Constructors*/
    //@{
    MultiSet()
    {
      size = 0;
      capacity = 0;
      lb_ = 0;
      ub_ = -1;

      values_ = NULL;
      index_ = NULL;
      occurrence_ = NULL;
    }

    MultiSet(const int lb, const int ub, const int* occ=NULL, const int* mocc=NULL)
    {
      if(!occ) {
	int tmp_occ[ub-lb+1];
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) {
	  tmp_occ[i] = 0;
	  tmp_mocc[i] = 4;
	}
	initialise(lb, ub, tmp_occ, tmp_mocc);
      } else if(!mocc) {
	int tmp_mocc[ub-lb+1];
	for(int i=0; i<ub-lb+1; ++i) 
	  tmp_mocc[i] = occ[i];
	initialise(lb, ub, occ, tmp_mocc);
      } else {
	initialise(lb, ub, occ, mocc);
      }
    }

    virtual ~MultiSet()
    {
      delete [] values_;
      occurrence_ += lb_;
      delete [] occurrence_;
      for(int v=lb_; v<=ub_; ++v)
	delete [] index_[v];
      index_ += lb_;
      delete [] index_;
    }

    void initialise(const int lb, const int ub, const int* occ, const int* mocc)
    {
      int span = ub-lb+1, i=0, j=0, v, o;

      lb_ = lb;
      ub_ = ub;

      index_ = new unsigned int*[span];
      occurrence_ = new unsigned int[span];
      capacity = 0;
      size = 0;

      for(v=0; v<span; ++v) {
	capacity += mocc[v];
	size += occ[v];
	occurrence_[v] = occ[v];
	index_[v] = new unsigned int[mocc[v]];
      }
      
      occurrence_ -= lb_;
      index_ -= lb_;
      values_ = new int[capacity];
      
      for(v=lb_; v<=ub_; ++v) {
	for(o=0; o<occ[v-lb_]; ++o) {
	  values_[i] = v;
	  index_[v][o] = i++;
	}
	for(o=occ[v-lb_]; o<mocc[v-lb_]; ++o) {
	  values_[size+j] = v;
	  index_[v][o] = size+j++;
	}
      }
    }
    //@}    


    inline void remove(const int elt)
    {
      --size;
      --occurrence_[elt];

      int v = values_[size];
      int i = index_[elt][occurrence_[elt]];

      index_[v][occurrence_[v]-1] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = values_[size];
      values_[size] = elt;
    }

    inline int pop()
    {
      int v = values_[--size];
      --occurrence_[v];
      return v;
    }

    inline int popHead()
    {
      const int elt = values_[0];

      --size;
      --occurrence_[elt];

      index_[values_[size]][occurrence_[values_[size]]-1] = 0;
      index_[elt][occurrence_[elt]] = size;

      values_[0] = values_[size];
      values_[size] = elt;

      return elt;
    }

    inline int head()
    {
      return *values_;
    }
    
    inline int back()
    {
      return values_[size-1];
    }

    inline void insert(const int elt)
    {
      int i = index_[elt][occurrence_[elt]];
      int v = values_[size]; 

      index_[v][occurrence_[v]] = i;
      index_[elt][occurrence_[elt]] = size;

      values_[i] = v;
      values_[size] = elt;

      ++occurrence_[elt];
      ++size;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
//     std::string getString() const {
//       std::string return_str = "{{";
//       bool first = true;
//       for(int v=lb_; v<=ub_; ++v) {
// 	if(occurrence_[v] > 0) {
// 	  if(first) first = false;
// 	  else return_str += ",";
// 	  return_str += ("("+toString((int)(occurrence_[v]))+"*"+toString(v)+")");
// 	}
//       }
//       return_str += "}}";
//       return return_str;
//     }

    std::ostream& display(std::ostream& os) const {
      os << "{{";
      bool first = true;
      for(int v=lb_; v<=ub_; ++v) {
	if(occurrence_[v] > 0) {
	  if(first) first = false;
	  else os << ",";
	  os << "(" << occurrence_[v] << "*" << v << ")" ;
	}
      }
      os << "}}";
      return os;
    }
    //@}
  };


  /**********************************************
   * IntStack
   **********************************************/
  /// Sparse set representation

  class IntStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    int *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    IntStack()
    {
      size = 0;
      capacity = 0;

      list_ = NULL;
      index_ = NULL;
    }

    IntStack(const int lb, const int ub, bool full=true)
    {
      initialise(lb, ub, full);
    }

    virtual ~IntStack()
    {
      delete [] list_;
      delete [] start_;
    }

    virtual void initialise(const int lb, const int ub, const bool full=true)
    {
      capacity = ub-lb+1;
      list_ = new int[capacity];
      start_ = new unsigned int[capacity];
      index_ = start_ - lb;      

      for(int i=lb; i<=ub; ++i) 
	{
	  index_[i] = i-lb;
	  list_[i-lb] = i;
	}
      
      size = (full ? capacity : 0);
    }

    void extend()
    {
      unsigned int new_capacity = capacity*2;
      int *aux_list = list_;
      list_ = new int[new_capacity];
      std::memcpy(list_, aux_list, capacity*sizeof(int));
      delete [] aux_list;

      unsigned int *aux_start = start_;
      start_ = new unsigned int[new_capacity];
      std::memcpy(start_, aux_start, capacity*sizeof(unsigned int));
      delete [] aux_start;

      capacity = new_capacity;
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool member(const int elt) const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline int next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline int operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline int& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void setTo(const int elt)
    {
      size=1;
      index_[*list_] = index_[elt];
      list_[index_[elt]] = *list_;
      *list_ = elt;
      index_[elt] = 0;
    }

    inline void erase(const int elt)
    {
      --size;
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
    }

    inline int next()
    {
      return list_[size];
    }

    inline int pop()
    {
      return list_[--size];
    }

    inline int popHead()
    {
      --size;
      index_[list_[size]] = 0;
      const int elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      return elt;
    }

    inline int head()
    {
      return *list_;
    }
    
    inline int back()
    {
      return list_[size-1];
    }

    inline void insert(const int elt)
    {
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];
      list_[size] = elt;
      index_[elt] = size;
      ++size;
    }

    inline void ordered_insert(const int elt)
    {
      // the first non-element goes where elt was
      index_[list_[size]] = index_[elt];
      list_[index_[elt]] = list_[size];

      int idx = size;
      while(idx && list_[idx-1] > elt) { // push every values greater than elt above elt
	list_[idx] = list_[idx-1];
	index_[list_[idx-1]] = idx;
	--idx;
      }

      list_[idx] = elt;
      index_[elt] = idx;
      ++size;
    }


    inline void revertTo(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i]] = i;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
//     std::string getString() const {
//       std::string return_str = "(";
//       if(size) return_str += toString(list_[0]);
//       for(unsigned int i=1; i<size; ++i)
// 	return_str += (" "+toString(list_[i]));
//       return_str += ")";
      
//       return return_str;
//     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };



  /**********************************************
   * Stack
   **********************************************/
  /// Sparse set representation

  template< class PTR_TYPE >
  class Stack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    PTR_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    Stack()
    {
      size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    Stack(Vector< PTR_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~Stack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void initialise(Vector< PTR_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back()->id - obj[0]->id + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new PTR_TYPE[capacity];
      offset = (obj.size ? obj[0]->id : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool member(const PTR_TYPE elt) const 
    {
      return index_[elt->id]<size;
    } 
    inline bool member(const int elt) const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline PTR_TYPE next(const PTR_TYPE elt) const
    {
      unsigned int idx = index_[elt->id]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline PTR_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline PTR_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline PTR_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void setTo(const PTR_TYPE elt)
    {
      int idx = elt->id;
      size=1;
      index_[(*list_)->id] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const PTR_TYPE elt)
    {
      int idx = elt->id;
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      PTR_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline PTR_TYPE next()
    {
      return list_[size];
    }

    inline PTR_TYPE pop()
    {
      return list_[--size];
    }

    inline PTR_TYPE popHead()
    {
      --size;
      index_[list_[size]->id] = 0;
      const PTR_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt->id] = size;
      return elt;
    }

    inline PTR_TYPE head()
    {
      return *list_;
    }
    
    inline PTR_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const PTR_TYPE elt)
    {
      int idx = elt->id;
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_insert(const PTR_TYPE elt)
    {
      int idx = elt->id;

      // the first non-element goes where elt was
      index_[list_[size]->id] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1]->id > elt->id) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1]->id] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revertTo(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i]->id] = list_[i]->id;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
//     std::string getString() const {
//       std::string return_str = "(";
//       if(size) return_str += toString(list_[0]);
//       for(unsigned int i=1; i<size; ++i)
// 	return_str += (" "+toString(list_[i]));
//       return_str += ")";
      
//       return return_str;
//     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * VarStack
   **********************************************/
  /// Sparse set representation

  template< class VAR_TYPE >
  class VarStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    VAR_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    /// values' indices
    unsigned int *index_;
    int offset;
    //unsigned int *start_;
    //@}

    /*!@name Constructors*/
    //@{
    VarStack()
    {
      size = 0;
      capacity = 0;
      list_ = NULL;
      offset = 0;
      index_ = NULL;
    }

    VarStack(Vector< VAR_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~VarStack()
    {
      delete [] list_;
      index_  += offset;
      delete [] index_;
    }

    void initialise(Vector< VAR_TYPE >& obj, const bool full=true)
    {
      assert((obj.size == 0) || ((unsigned int)(obj.back().id() - obj[0].id() + 1) == obj.size));

      capacity = (obj.size ? obj.size : obj.capacity);
      list_ = new VAR_TYPE[capacity];
      offset = (obj.size ? obj[0].id() : 0);
      index_ = new unsigned int[capacity];
      index_ -= offset;
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  index_[i+offset] = i;
	  list_[i] = obj[i];
	}
      
      size = (full ? capacity : 0);
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool member(const VAR_TYPE elt) const 
    {
      return index_[elt.id()]<size;
    } 
    inline bool member(const int elt) const 
    {
      return index_[elt]<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline VAR_TYPE next(const VAR_TYPE elt) const
    {
      unsigned int idx = index_[elt.id()]+1;
      return (idx < size ? list_[idx] : elt);
    }
    inline VAR_TYPE next(const int elt) const
    {
      unsigned int idx = index_[elt]+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline VAR_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline VAR_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void setTo(const VAR_TYPE elt)
    {
      int idx = elt.id();
      size=1;
      index_[(*list_).id()] = index_[idx];
      list_[index_[idx]] = *list_;
      *list_ = elt;
      index_[idx] = 0;
    }

    inline void remove(const VAR_TYPE elt)
    {
      int idx = elt.id();
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline void remove(const int idx)
    {
      VAR_TYPE elt = list_[index_[idx]];
      --size;
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
    }

    inline VAR_TYPE next()
    {
      return list_[size];
    }

    inline VAR_TYPE pop()
    {
      return list_[--size];
    }

    inline VAR_TYPE popHead()
    {
      --size;
      index_[list_[size].id()] = 0;
      const VAR_TYPE elt = *list_;
      *list_ = list_[size];
      list_[size] = elt;
      index_[elt.id()] = size;
      return elt;
    }

    inline VAR_TYPE head()
    {
      return *list_;
    }
    
    inline VAR_TYPE back()
    {
      return list_[size-1];
    }

    inline void add(const VAR_TYPE elt)
    {
      int idx = elt.id();
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];
      list_[size] = elt;
      index_[idx] = size;
      ++size;
    }

    inline void ordered_insert(const VAR_TYPE elt)
    {
      int idx = elt.id();

      // the first non-element goes where elt was
      index_[list_[size].id()] = index_[idx];
      list_[index_[idx]] = list_[size];

      int rank = size;
      while(idx && list_[rank-1].id() > elt.id()) { // push every values greater than elt above elt
	list_[rank] = list_[rank-1];
	index_[list_[rank-1].id()] = rank;
	--rank;
      }

      list_[rank] = elt;
      index_[idx] = rank;
      ++size;
    }

    inline void revertTo(const int level)
    {
      size = level;
    }

    inline void index()
    {
      for(unsigned int i=0; i<capacity; ++i)
	index_[list_[i].id()] = list_[i].id();
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
//     std::string getString() const {
//       std::string return_str = "(";
//       if(size) return_str += toString(list_[0]);
//       for(unsigned int i=1; i<size; ++i)
// 	return_str += (" "+toString(list_[i]));
//       return_str += ")";
      
//       return return_str;
//     }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * IdStack
   **********************************************/
  /// Sparse set representation
  template< class DATA_TYPE >
  class IdStack 
  {
  public:

    /*!@name Parameters*/
    //@{
    /// list of values
    DATA_TYPE *list_;
    /// current max capacity
    unsigned int capacity;
    /// current size
    unsigned int size;
    //@}

    /*!@name Constructors*/
    //@{
    IdStack()
    {
      size = 0;
      capacity = 0;
      list_ = NULL;
    }

    IdStack(Vector< DATA_TYPE >& obj, bool full=true)
    {
      initialise(obj, full);
    }

    virtual ~IdStack()
    {
      delete [] list_;
    }

    void initialise(Vector< DATA_TYPE >& obj, const bool full=true)
    {
      capacity = obj.size;
      size = (full ? capacity : 0);
      list_ = new DATA_TYPE[capacity];
      for(unsigned int i=0; i<capacity; ++i) 
	{
	  list_[i] = obj[i];
	  list_[i].set_stack_id(i);
	}
    }
    //@}    

    /*!@name Accessors*/
    //@{  
    inline bool member(const DATA_TYPE elt) const 
    {
      return elt.get_stack_id()<size;
    } 
  
    inline bool empty()const 
    {
      return !size;
    } 

    inline DATA_TYPE next(const DATA_TYPE elt) const
    {
      unsigned int idx = elt.get_stack_id()+1;
      return (idx < size ? list_[idx] : elt);
    }
    
    inline DATA_TYPE operator[](const unsigned int idx) const
    {
      return list_[idx];
    }

    inline DATA_TYPE& operator[](const unsigned int idx)
    {
      return list_[idx];
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline void fill()
    {
      size = capacity;
    }

    inline void clear()
    {
      size = 0;
    }
  
    inline void setTo(DATA_TYPE elt)
    {
      int idx = elt.get_stack_id();
      size=1;
      
      (*list_).set_stack_id(idx);
      elt.set_stack_id(0);

      list_[idx] = *list_;
      *list_ = elt;
    }

    inline void remove(DATA_TYPE elt)
    {
      int idx = elt.get_stack_id();

      list_[--size].set_stack_id(idx);
      elt.set_stack_id(size);
      
      list_[idx] = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE pop()
    {
      return list_[--size];
    }

    inline DATA_TYPE popHead()
    {
      --size;
      const DATA_TYPE elt = *list_;

      list_[size].set_stack_id(0);
      elt.set_stack_id(size);

      *list_ = list_[size];
      list_[size] = elt;
    }

    inline DATA_TYPE head()
    {
      return *list_;
    }
    
    inline DATA_TYPE back(const int offset=0)
    {
      return list_[size-1+offset];
    }

    inline void add(DATA_TYPE elt)
    {
      int idx = elt.get_stack_id();
      
      elt.set_stack_id(size);
      list_[size].set_stack_id(idx);

      list_[idx] = list_[size];
      list_[size] = elt;
      
      ++size;
    }

    inline void revertTo(const int level)
    {
      size = level;
    }

    std::ostream& display(std::ostream& os) const {
      os << "(";
      if(size) os << list_[0];
      for(unsigned int i=1; i<size; ++i)
	os << " " << list_[i];
      os << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * Node
   **********************************************/
  /// Node of Multilist

  template < class DATA_TYPE >
  class Node {
  public:

    /*!@name Parameters*/
    //@{
    int prev;
    int next;
    DATA_TYPE elt;
    //@}

    Node(DATA_TYPE e=(DATA_TYPE(0)), int p=-1, int n=-1) {elt = e; prev=p; next=n;}
    inline operator const DATA_TYPE() const { return elt; }

    //inline Node& operator=(const Node& nd) const { prev = nd.prev; return elt; }

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
//     std::string getString() const {
//       std::string return_str = "("+toString(prev)+"|"+toString(elt)+"|"+toString(next)+")";
//       return return_str;
//     }
    std::ostream& display(std::ostream& os) const {
      os << "(" << prev << "|" << elt << "|" << next << ")";
      return os;
    }
    //@}
  };


  /**********************************************
   * Multilist
   **********************************************/
  /// List with multiple entry points

  template < class DATA_TYPE, int NUM_HEAD >
  class MultiList 
  {

  public:
    /*!@name Parameters*/
    //@{
    Vector< Node< DATA_TYPE > > data;
    int head[NUM_HEAD+1]; // index of the first element of the list
    unsigned int degree; // number of element in the entire list
    //@}

    /*!@name Constructors*/
    //@{
    MultiList() {
      data.initialise(0, std::max(2*NUM_HEAD,16));
      for(int k=0; k<=NUM_HEAD; ++k) {
	head[k] = k; 
	Node< DATA_TYPE > x((DATA_TYPE)0, k-1, k+1); 
	data.add(x);
      }
      degree=0;
    }
    //@}

    /*!@name List Manipulation*/
    //@{
    inline int create(DATA_TYPE elt, const int k=0) {
      Node<DATA_TYPE> x(elt);
      data.add(x);
      insert(data.size-1, k);
      return data.size-1;
    }

    inline void insert(const int idx, const int k=0) {
      int succ = data.stack_[head[k]].next;
      data.stack_[succ].prev = idx;
      data.stack_[head[k]].next = idx;
      data.stack_[idx].next = succ;
      data.stack_[idx].prev = head[k];
      ++degree;
    }

    inline void erase(const int idx, const int k=0) {
      int succ = data.stack_[idx].next;
      int prec = data.stack_[idx].prev;      
      data.stack_[succ].prev = prec;
      data.stack_[prec].next = succ;
      --degree;
    }

    inline Node<DATA_TYPE>& first(const int h) {
      //int idx = ;
      return data.stack_[head[h]];
    }

    inline Node<DATA_TYPE>& last(const int h) {
      //int idx = data.stack_[head[h+1]].prev;
      return data.stack_[data.stack_[head[h+1]].prev];
    }

    inline DATA_TYPE pop(const int h) {
      int idx = data.stack_[head[h+1]].prev;
      erase(idx, h);
      return data.stack_[idx].elt;
    }

    inline bool empty(const int h) {
      return (data.stack_[head[h]].next > NUM_HEAD);
    }

    inline bool next(Node< DATA_TYPE >& node) const {
      unsigned int idx_next;
      do { 
	idx_next = node.next;
	if(idx_next == NUM_HEAD) return false;
	node = data.stack_[idx_next];
      } while(!node.elt);
      return true;
    }

    inline bool prev(Node< DATA_TYPE >& node) const {
      unsigned int idx_prev;
      do { 
	idx_prev = node.prev;
	if(!idx_prev) return false;
	node = data.stack_[idx_prev];
      } while(!node.elt);
      return true;
    }
    //@}

    /*!@name Miscellaneous*/
    //@{
    /// Print on out stream
//     std::string getString() const {
//       std::string return_str = "";
//       for(unsigned int i=0; i<data.size; ++i)
// 	return_str += toString(data[i]);
//       return return_str;
//     }
    std::ostream& display(std::ostream& os) const {
      for(unsigned int i=0; i<data.size; ++i)
	os << data[i];
      return os;
    }

    void debug_print(std::ostream& o) const 
    {
      for(unsigned int i=head[0]; i<data.size; ++i)
	data.stack_[i].print(o);
    }
    void flat_print(std::ostream& o) const 
    {
      o << "[ ";
      for(unsigned int i=NUM_HEAD+1; i<data.size; ++i)
	o << data.stack_[i].elt << " ";
      o << "] ";
    }
    void print(std::ostream& o) const 
    {
      for(int k=0; k<NUM_HEAD; ++k) {
	o << "[ ";

	Node< DATA_TYPE > nd = data.stack_[head[k]];
	while(next(nd)) {
	  o << (DATA_TYPE)nd << " ";
	}
	o << "] ";
      }
    }
    //@}

  };


  /**********************************************
   * BitSet
   **********************************************/

  /*! \class BitSet
    \brief A representation of sets using a vector of bits.

    The sets have a static capacity. 
    Template are used so that ReversibleWords can used instead of unsigned int
  */
  template< class WORD_TYPE, class FLOAT_TYPE >
  class Bitset
  { 

  public:
    
    /*!@name Class attributes*/
    //@{
    static const WORD_TYPE empt = 0;
    static const WORD_TYPE full = ~0;
    static const unsigned int EXP = (sizeof(empt) == 4 ? 5 /*32 bits*/ : 6 /*64 bits*/);
    static const unsigned int size_word_bit = (1 << EXP);
    static const unsigned int size_word_byte = (size_word_bit >> 3);
    static const unsigned int CACHE = (size_word_bit - 1);
    static const unsigned int LASTCHAR = (size_word_bit - 8);
    static const unsigned int mantissa = (sizeof(empt) == 4 ? 23 /*32 bits*/ : 52 /*64 bits*/);
    static const unsigned int float_offset = (sizeof(empt) == 4 ? (0x7f) /*32 bits*/ : (0x3ff) /*64 bits*/);
    static const WORD_TYPE mask_first_char = 0xff;
    static const WORD_TYPE mask_last_char = (mask_first_char << ((size_word_bit) - 8));
    //@}

    /*!@name Parameters*/
    //@{
    /// index of the first word used to represent the set
    int pos_words;
    /// 1 + index of the last word used to represent the set
    int neg_words;
    /// A vector of bits 
    WORD_TYPE* table;
    //@}

    Bitset()
    {
      initialise();
    }

    void initialise()
    {
      pos_words = 0;
      neg_words = 0;
      table = NULL;
    }

    Bitset(const int sz, const int* elt) 
    {
      int lb = INT_MAX;
      int ub = INT_MIN;
      for(int i=0; i<sz; ++i) {
	if(elt[i] > ub) ub = elt[i];
	if(elt[i] < lb) lb = elt[i];
      }

      initialise(lb,ub,empt);

      for(int i=0; i<sz; ++i) 
	insert( elt[i] );
    }

    Bitset(const int lb, const int ub, const WORD_TYPE p)//, WORD_TYPE *pool=NULL) 
    {
      initialise(lb,ub,p,NULL);
    }

    void reinitialise(const int lb, const int ub, const WORD_TYPE p) 
    {
      table += neg_words;
      delete [] table;
      initialise(lb, ub, p, NULL);
    }

    void initialise(const int sz, const WORD_TYPE p) 
    {
      pos_words = sz;
      neg_words = 0;

      if( sz>=0 ) {
	table = new WORD_TYPE[pos_words];
	for(unsigned int i=0; i<pos_words; ++i) 
	  table[i]=p;
      } else table = NULL;
    }

    void initialise(const int lb, const int ub, const WORD_TYPE p, WORD_TYPE *pool=NULL) 
    {
      neg_words = (lb >> EXP);
      pos_words = (ub >> EXP)+1;
      if(pool==NULL) table = new WORD_TYPE[pos_words-neg_words];
      else table = pool;
      for(int i=0; i<pos_words-neg_words; ++i) 
	table[i]=p;
      table[pos_words-neg_words-1] &= 
	(p >> (size_word_bit-1-(ub & CACHE)));
      table[0] &= (p << (lb & CACHE));
      table -= neg_words;
    }

    void initialise(Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      pos_words = s.pos_words;
      neg_words = s.neg_words;
  
      table = new WORD_TYPE[pos_words-neg_words];
      table -= neg_words;
      for(unsigned int i=neg_words; i<pos_words; ++i) 
	table[i].initialise(s.table+i, s.size(i));
    }

    Bitset(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      clone( s );
    }

    void clone(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = new WORD_TYPE[pos_words-neg_words];
      memcpy(table, s.table+neg_words,
	     size_word_byte*(pos_words-neg_words));
      table -= neg_words;
    }

    void pointTo(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      neg_words = s.neg_words;
      pos_words = s.pos_words;
      table = s.table;
    }

    void pointTo(WORD_TYPE *t)
    {
      neg_words = 0;
      pos_words = 1;
      table = t;
    }

    void copy(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int k, j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      for( k=neg_words; k<j; ++k)
	table[k] = empt;
      for( k=i; k<pos_words; ++k)
	table[k] = empt;
      if( i>j )
	memcpy(table+j,s.table+j,size_word_byte*(i-j));
    }

    virtual ~Bitset() 
    {
      table += neg_words;
      delete [] table; 
    }

    void destroy() 
    {
      table += neg_words;
      neg_words = 0;
      delete [] table; 
      table = NULL;
    }

    bool isBuilt()
    {
      return (table != NULL);
    }

    inline void swap(Bitset<WORD_TYPE,FLOAT_TYPE>& s)
    {
      WORD_TYPE *aux = s.table;
      s.table = table;
      table = aux;
    }

    inline int minimum_element(int idx, WORD_TYPE v) const
    {
      while(v == 0) {
	if( ++idx >= pos_words )
	  return NOVAL;
	v = table[idx];
      }

      union {FLOAT_TYPE f; WORD_TYPE i; } t;
      WORD_TYPE b = v & -v;

      t.f = (FLOAT_TYPE)b; // cast the least significant bit in v to a float
      b = t.i >> mantissa;
      return b + idx * size_word_bit - float_offset;
    }

    /*!
      Minimum element in the set [O(N/8)]. 
    */
    inline int min() const
    { 
      int idx = neg_words;
      WORD_TYPE v = table[idx];
      return minimum_element(idx,v);
    }

    /*!
      Maximum element in the set [O(N/8)]. 
    */
    inline int max() const
    { 
      WORD_TYPE tab;
      int i=pos_words, j, k;
    
      while( i-- > neg_words )
	if( (tab = table[i]) ) {
	  j = size_word_byte;
	  while( j-- ) {
	    if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 ) 
	      return ( (i<<EXP)+(j<<3)+k );	
	    tab = (tab << 8);
	  }
	}
      return NOVAL;
    }

    inline void  erase(const int elt) 
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline void  fastErase(const int elt) 
    {
      table[(elt >> EXP)] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline void  wordErase(const int elt) 
    {
      table[neg_words] &= (full ^ ((WORD_TYPE)1 << (elt & CACHE)));
    }

    inline int next(const int elt) const {
      int idx = ((elt+1) >> EXP);
      WORD_TYPE v = (table[idx] & (full << ((elt+1) & CACHE)));
      return minimum_element(idx,v);
    }

    inline int prev(const int elt) const {
      WORD_TYPE tab;
      int i = ((elt-1) >> EXP);

      if( i >= neg_words ) {
	int e = ((elt-1) & CACHE), k;
	int j = 1+(e >> 3);

	if( (tab = ((table[i] & (full >> (CACHE - e))) << ((4-j) << 3))) ) 
	  while( j-- ) {
	    if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
	      return ( (i<<EXP)+(j<<3)+k );
	    tab = (tab << 8);
	  }
	while( i-- > neg_words ) 
	  if( (tab = table[i]) ) {
	    j = size_word_byte;
	    while( j-- ) {
	      if( (k = getlast[(tab & mask_last_char) >> LASTCHAR]) >= 0 )
		return ( (i<<EXP)+(j<<3)+k );
	      tab = (tab << 8);
	    }
	  }
      }
      return NOVAL;
    }

    inline void xorTo(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	s.table[i] ^= table[i];
    }

    inline void fastXorTo(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	s.table[i] ^= table[i];
    }

    inline void xorWith(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] ^= s.table[i];
    }

    inline void fastXorWith(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = pos_words;
      while( i-- > neg_words )
	table[i] ^= s.table[i];
    }

    inline void unionWith(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] |= s.table[i];
    }

    inline void unionWith(const int s) 
    {
      if(pos_words>0 && neg_words<=0) table[0] |= s;
    }

    inline void unionTo(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      s.unionWith( *this );
    }

    inline void intersectWith(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	table[k] = empt;
      }
      while( k > j ) {
	--k;
	table[k] &= s.table[k];
      }
      while( k-- > neg_words )
	table[k] = empt;
    }

    inline void intersectWith(const int s) 
    {
      int i=pos_words;
      while(i-- > 1) table[i] = empt;
      i = 0;
      while(i-- > neg_words) table[i] = empt;
      if(pos_words>0 && neg_words<=0) table[0] &= s;
    }

    inline void intersectTo(Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.intersectWith( *this );
    }

    inline void setminusWith (const Bitset<WORD_TYPE,FLOAT_TYPE>& s) 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	table[i] &= (~(s.table[i]));
    }

    inline void setminusTo (Bitset<WORD_TYPE,FLOAT_TYPE>& s) const
    {
      s.setminusWith( *this );
    }

    inline void xorWith(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      xorWith(*s);
    }

    inline void unionWith  (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      unionWith(*s);
    }

    inline void intersectWith(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      intersectWith(*s);
    }

    inline void setminusWith (const Bitset<WORD_TYPE,FLOAT_TYPE>* s) 
    {
      setminusWith(*s);
    }
    inline void unionTo  (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->unionWith(*this);
    }

    inline void intersectTo(Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->intersectWith(*this);
    }

    inline void setminusTo (Bitset<WORD_TYPE,FLOAT_TYPE>* s) const
    {
      s->setminusWith(*this);
    }

    inline bool equal(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      int i=pos_words;
      int j=s.pos_words;
      int k;
  
      while( j > i )
	if( s.table[--j] ) return false; 
      while( i > j )
	if( table[--i] ) return false;

      j=neg_words;
      k=s.neg_words;

      while( j > k )
	if( s.table[k++] ) return false;
      while( k > j )
	if( table[j++] ) return false;

      while( i-- > j )
	if( table[i] != s.table[i] ) return false;
  
      return true;
    }
    
    inline bool includes(const WORD_TYPE s) const 
    {
      return( pos_words && neg_words<1 && (table[0] & s) == s );
    }

    inline bool included(const WORD_TYPE s) const 
    {
      bool inc = true;
      int k = pos_words;
      if(neg_words>0 || pos_words<1) {
	while(k>neg_words && inc) inc = !(table[--k]);
      } else {
	while(k>1 && inc) inc = !(table[--k]);
	inc = ((table[--k] & s) == table[0]);
	while(k>neg_words && inc) inc = !(table[--k]);
      }
      return inc;
    }

    inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	if( table[k] ) return false;
      }
      while( k > j ) {
	--k;
	if( table[k] != (table[k] & s.table[k]) ) return false;
      }
      while( k-- > neg_words )
	if( table[k] ) return false;
      return true;
    }

    inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      int k = pos_words;
      while( k > i ) {
	--k;
	if( s.table[k] ) return false;
      }
      while( k > j ) {
	--k;
	if( s.table[k] != (table[k] & s.table[k]) ) return false;
      }
      while( k-- > neg_words )
	if( s.table[k] ) return false;
      return true;
    }

    inline bool included(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return included( *s );
    }

    inline bool includes(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return includes( *s );
    }

    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>* s) const 
    {
      return intersect( *s );
    }

    inline bool included(const int lb, const int ub) const 
    {
      int neg_int = lb >> EXP;
      int pos_int = ub >> EXP;
      int k = pos_words;
      while( k > pos_int ) 
	if( table[--k] ) return false;
      k = neg_words;
      while( k < neg_int )
	if( table[k++] ) return false;
      if(neg_int == pos_int) {
	k = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
	return (k & table[neg_int]) == table[neg_int];
      } else {
	return ((((full << (lb & CACHE)) & table[neg_int]) == table[neg_int]) &&
		(((full >> (CACHE - (ub & CACHE))) & table[pos_int]) == table[pos_int]));
      }
    }

    inline bool includes(const int lb, const int ub) const
    {
      int neg_int = lb >> EXP;
      int pos_int = ub >> EXP;
      int k = pos_int-1;
      unsigned int u, l;
      while( k > neg_int ) {
	if( table[k] != full ) return false;
	--k;
      }
      if(neg_int == pos_int) {
	u = ((full << (lb & CACHE)) & (full >> (CACHE - (ub & CACHE))));
	return (u & table[neg_int]) == u;
      } else {
	u = (full >> (CACHE - (ub & CACHE)));
	l = (full << (lb & CACHE));
	return (((l & table[neg_int]) == l) &&
		((u & table[pos_int]) == u));
      }
    }

    /*!
     * Returns the number of bits set in v.
     * For a derivation of this algorithm, see
     * "Algorithms and data structures with applications to 
     *  graphics and geometry", by Jurg Nievergelt and Klaus Hinrichs,
     *  Prentice Hall, 1993.
     */
    inline unsigned int word_size(WORD_TYPE v) const 
    {
      v = v - ((v >> 1) & (WORD_TYPE)~(WORD_TYPE)0/3);                           // temp
      v = (v & (WORD_TYPE)~(WORD_TYPE)0/15*3) + ((v >> 2) & (WORD_TYPE)~(WORD_TYPE)0/15*3);      // temp
      v = (v + (v >> 4)) & (WORD_TYPE)~(WORD_TYPE)0/255*15;                      // temp
      return (WORD_TYPE)(v * ((WORD_TYPE)~(WORD_TYPE)0/255)) >> (sizeof(v) - 1) * CHAR_BIT; // count
    }
	
    inline unsigned int size() const 
    {  
      int i=pos_words;
      unsigned int c=0;
      WORD_TYPE v;
      while( i-- > neg_words ) 
	if( (v = table[i]) ) 
	  c += word_size(v);
      return c;  
    }

    inline unsigned int wordSize() const 
    {  
      unsigned int v, c=0;
      if( (v = table[neg_words]) ) 
	c = word_size(v);
      return c;  
    }

    inline unsigned int size( const int i ) const
    {  
      unsigned int v, c=0;
      if( (v = table[i]) ) 
	c = word_size(v);
      return c;  
    }


    /*!
      Check if element elt belong to the set [O(1)]
    */
    inline  bool member(const int elt)const 
    {
      int i = (elt >> EXP);
      return ( (i >= neg_words) && 
	       (i <  pos_words) && 
	       (table[i] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool fastMember(const int elt)const 
    {
      return ( (table[(elt >> EXP)] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }

    inline  bool wordMember(const int elt)const 
    {
      return ( (table[neg_words] & ((WORD_TYPE)1 << (elt & CACHE))) );
    }
    /*!
      Insert element elt into the set [O(1)]
    */

    inline  void insert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) ) 
	table[i] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fastInsert(const int elt)
    {
      table[(elt >> EXP)] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline   void wordInsert(const int elt)
    {
      table[neg_words] |= ((WORD_TYPE)1 << (elt & CACHE));
    }

    /*!
      Insert element elt into the set or erase it if it is already member [O(1)]
    */
    inline  void invert(const int elt)
    {
      int i = (elt >> EXP);
      if( (i >= neg_words) && 
	  (i <  pos_words) )
	table[i] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void fastInvert(const int elt)
    {
      table[(elt >> EXP)] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    inline  void wordInvert(const int elt)
    {
      table[neg_words] ^= ((WORD_TYPE)1 << (elt & CACHE));
    }

    /*!
      Return true iff the set is empty [O(N/32)]
    */
    inline  bool empty()const
    { 
      int i = pos_words;
      while( i-- > neg_words ) 
	if(table[i]) return false;
      return true;  
    }

    /*!
      Return true iff the calling object intersect s [O(N/32)]
    */
    inline bool intersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s)const
    {
      int i = (pos_words > s.pos_words ? s.pos_words : pos_words);
      int j = (neg_words < s.neg_words ? s.neg_words : neg_words);
      while( i-- > j )
	if(table[i] & s.table[i]) return true;
      return false;
    }

    /*!
      Return true iff the calling object intersect s (s is assumed to be a bitset in {0,..,31}) [O(N/32)]
    */
    inline bool intersect(const int s) const
    {
      return(pos_words && neg_words<1 && (table[0]&s));
    }

    inline bool wordIntersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s) const 
    {
      return ( table[neg_words] & s.table[neg_words] ) ;
    }

    inline  bool fastIntersect(const Bitset<WORD_TYPE,FLOAT_TYPE>& s, int& idx)const
    {
      if( table[idx] & s.table[idx] ) return true; 
      if( pos_words > neg_words ) {
	idx = pos_words;
	while( idx > neg_words ) {
	  --idx;
	  if(table[idx] & s.table[idx]) return true;
	}
      }
      return false;
    }

    /*!
      Return true iff the calling object intersect [lo..up] [O(N/32)]
    */
    inline   bool intersect(const int lb, const int ub)const
    {
      int i = (ub >> EXP);
      int j = (lb >> EXP);

      if( i < neg_words || j >= pos_words ) 
	return false;

      WORD_TYPE masked_lb = (full << (lb & CACHE));
      WORD_TYPE masked_ub = (full >> (CACHE - (ub & CACHE)));

      if( i == j ) 
	if( table[i] & (masked_lb & masked_ub) ) return true;
	else return false;

      if( i >= pos_words )
	i = pos_words-1;
      else if( table[i--] & masked_ub ) return true;
  
      if( j < neg_words ) 
	j = neg_words;
      else if( table[j++] & masked_lb ) return true;

      while( i >= j )
	if(table[i--]) return true;
      return false;
    }

    /*!
      Increment by x all elements in the set.
      Any element greater or equal than the capacity 
      of the set minus x is erased [O(N/32)]
    */
    inline void increment(const int v)
    {
      int step = (v >> EXP); 
      int i = pos_words;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = neg_words+step;
      WORD_TYPE mask = ((WORD_TYPE)~0 << f);
      while( --i > j ) 
	table[i] = ((table[i-step] << e) | ((table[i-step-1] & mask) >> f));
      if( i >= neg_words+step ) table[i] = (table[i-step] << e);
      while( i > neg_words ) table[--i] = 0;
    }

    /*!
      Decrement by x all elements in the set.
      Any element lower than x is erased [O(N/32)]
    */
    inline  void decrement(const int v)
    {
      int step = (v >> EXP); 
      int i = neg_words-1;
      int e = (v & CACHE);
      int f = size_word_bit-e;
      int j = pos_words-step-1;
      WORD_TYPE mask = ((WORD_TYPE)~0 >> e);
      while( ++i < j )
	table[i] = ((table[i+step] >> e) | ((table[i+step+1] & mask) << f));
      if( i < pos_words-step ) table[i] = (table[i+step] >> e);
      while( ++i < pos_words ) table[i] = 0;
    }

    /*!
      Changes every value to its arythmetic negation
    */
    inline void negate( Bitset<WORD_TYPE,FLOAT_TYPE>& s )
    {
      int i = (pos_words > -s.neg_words ? -s.neg_words : pos_words);
      int j = (neg_words < -s.pos_words ? -s.pos_words : neg_words);

      unsigned int a; 
      WORD_TYPE mask, v, aux, rest = ( i < pos_words && (table[i] & 1) );

      while( i-- > j ) {
	aux = (table[i] & 1);
	v = (table[i] >> 1);
	mask = ~0;         
	a = sizeof(v) * CHAR_BIT; // bit size; must be power of 2 
	while ((a >>= 1) > 0) 
	  {
	    mask ^= (mask << a);
	    v = ((v >> a) & mask) | ((v << a) & ~mask);
	  }	
	s.table[-i-1] = (v | rest);
	rest = aux;
      }
      if(rest)
	s.table[i+1] |= rest;
    }

    /*!
      Insert all elements between 0 to capacity [O(N/32)]
    */
    inline  void fill()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = full;
    }

    /*!
      Erase all elements [O(N/32)]
    */
    inline  void clear()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] = empt;
    }

    /*!
      Erase all elements but v [O(N/32)]
    */
    inline  void setTo( const int v )
    {
      int i, j = (v >> EXP);
      for(i=neg_words; i<j; ++i)
	table[i] = empt;
      table[j] = ((WORD_TYPE)1 << v);
      for(i=j+1; i<pos_words; ++i) 
	table[i] = empt;
    }

    /*!
      flip all elements [O(N/32)]
    */
    inline  void flip()
    {
      int i = pos_words;
      while( i > neg_words )
	table[--i] ^= full;
    }

    /*!
      Erase all elements strictly lower than l [O(N/32)]
    */
    inline void setMin(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word >= neg_words )
	if( ith_word <  pos_words ) {      
	  int i=ith_word;
	  while( i-- > neg_words ) table[i]=0;
	  table[ith_word] &= (full << (bound & CACHE));
	} else clear();
    }

    /*!
      Erase all elements strictly greater than u [O(N/32)]
    */
    inline void setMax(const int bound)
    {
      int ith_word=(bound >> EXP);
      if( ith_word <  pos_words )
	if( ith_word >= neg_words ) {
	  int i=pos_words;
	  while( --i > ith_word ) table[i]=0;
	  table[ith_word] &= (full >> (CACHE - (bound & CACHE)));
	} else clear();
    }

    /*!
      Erase all elements in the interval [l..u] [O(N/32)]
    */
    inline  void removeInterval(const int lb, const int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = 0;
	WORD_TYPE masked_ub = 0;

	if( lb_word >= neg_words ) 
	  masked_lb = (full >> (CACHE - (lb & CACHE) + 1));
	if( ub_word < pos_words ) 
	  masked_ub = ((full-1) << (ub & CACHE) );

	if( lb_word == ub_word ) 
	  table[lb_word] &= (masked_lb | masked_ub);
	else {
	  table[lb_word] &= masked_lb;
	  table[ub_word] &= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = 0;
	}
      }
    }

    /*!
      Add all elements in the interval [l..u] [O(N/32)]
    */
    inline void addInterval(int lb, int ub)
    {
      if( lb <= ub ) {
	int lb_word = lb >> EXP;
	int ub_word = ub >> EXP;
    
	lb_word = ( lb_word < neg_words ? neg_words : lb_word );
	ub_word = ( ub_word >= pos_words ? pos_words-1 : ub_word );

	WORD_TYPE masked_lb = full;
	WORD_TYPE masked_ub = full;
	if( lb_word >= neg_words ) 
	  masked_lb ^= (full >> (CACHE - (lb & CACHE) + 1));
	if( ub_word < pos_words ) 
	  masked_ub ^= ((full-1) << (ub & CACHE) );
    
	if( lb_word == ub_word ) 
	  table[lb_word] |= (masked_lb & masked_ub);
	else {
	  table[lb_word] |= masked_lb;
	  table[ub_word] |= masked_ub;
	  while( --ub_word > lb_word )
	    table[ub_word] = full;
	}
      }
    }

    inline bool operator[](const int i)
    {
      return fastMember(i);
    }

//     std::string getString() const {  
//       std::string return_str = "{";
//       if( !empty() ) {
// 	int last = NOVAL, cur=min(), aft;

// 	bool flag=false;
// 	do{
// 	  aft = next(cur);

// 	  if(aft != cur+1 || cur != last+1) {
// 	    if( flag )
// 	      return_str += ",";
// 	    return_str += toString(cur);
// 	    flag = true;
// 	  } else if(flag) {
// 	    return_str += "..";
// 	    flag = false;
// 	  }
// 	  last = cur;
// 	  cur = aft;
// 	} while( cur != NOVAL && cur != last );
//       }
//       return_str += "}";
//       return return_str;
//     }

    std::ostream& display(std::ostream& os) const {
      os << "{";
      if( !empty() ) {
	int last = NOVAL, cur=min(), aft;

	bool flag=false;
	do{
	  aft = next(cur);

	  if(aft != cur+1 || cur != last+1) {
	    if( flag )
	      os << ",";
	    os << cur;
	    flag = true;
	  } else if(flag) {
	    os << "..";
	    flag = false;
	  }
	  last = cur;
	  cur = aft;
	} while( cur != NOVAL && cur != last );
      }
      os << "}";
      return os;
    }

//     void  print(std::ostream & os) const 
//     {
//       //simplePrint();
//       if( table )
// 	print( os, "{,}" );
//     }

//     void  print(std::ostream & os, 
// 		const char *delim) const
//     {
//       os << delim[0];
//       if( !empty() ) {
// 	int last = NOVAL, cur=min(), aft;

// 	bool flag=false;
// 	do{
// 	  aft = next(cur);

// 	  if(aft != cur+1 || cur != last+1) {
// 	    if( flag )
// 	      os << delim[1];
// 	    os << cur;
// 	    flag = true;
// 	  } else if(flag) {
// 	    os << "..";
// 	    flag = false;
// 	  }
// 	  last = cur;
// 	  cur = aft;
// 	} while( cur != NOVAL && cur != last );
//       }
//       os << delim[2];
//     }

    void  printBits(std::ostream & os) const 
    {
      os << "[";
      for(int i=neg_words; i<pos_words; ++i)
 	showUint( table[i], os );
      os << "]";
    }

  };


  /**********************************************
   * InlinedBitSet
   **********************************************/
  /// The data array is stored directly in the object

  template< class WORD_TYPE, class FLOAT_TYPE >
  class InlinedBitset : public Bitset<WORD_TYPE,FLOAT_TYPE> 
  {

  public:

    /**@name Parameters*/
    //@{    
    WORD_TYPE data[0];
    //@}
    
    static InlinedBitset<WORD_TYPE,FLOAT_TYPE>* make_new(const int lb, const int ub, const WORD_TYPE p)
    {
      int n_w = (lb >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP);
      int p_w = (ub >> Bitset<WORD_TYPE,FLOAT_TYPE>::EXP)+1;
      int size = p_w - n_w;
      
      void* mem = malloc(sizeof(InlinedBitset<WORD_TYPE,FLOAT_TYPE>) + sizeof(WORD_TYPE)*(size));
      return new (mem) InlinedBitset<WORD_TYPE,FLOAT_TYPE>(lb, ub, p);
    }

    InlinedBitset(const int lb, const int ub, const WORD_TYPE p)
    {
      initialise(lb,ub,p,&(data[0]));
    }
    
  };

  typedef Bitset< unsigned long long int, double > Bitset64;
  typedef Bitset< unsigned int, float > Bitset32;
  typedef InlinedBitset< unsigned long long int, double > iBitset64;
  typedef InlinedBitset< unsigned int, float > iBitset32;

#ifdef _BIT64
  
  typedef Bitset64 BitSet;

#else

  typedef Bitset32 BitSet;

#endif


  /**********************************************
   * Queue
   **********************************************/
  /// List of constraints for computing the GAC closure (fifo).
  
  class Queue {
  public:
    
    /**@name Parameters*/
    //@{
    //BitSet _set_;
    int *next;
    int _head;
    int _tail;
    //@}

    /**@name Constructors*/
    //@{
    Queue() {
      next = NULL;
    }
    ~Queue()
    {
      delete [] next;
    }
    void initialise(const int n)
    {
      next = new int[n+1];
      std::fill( next, next+n, NOVAL );
      _head = _tail = n;
      next[n] = n;
    }
    //@}
  
    /*!@name Accessors*/
    //@{
    /// Add an element to the queue
    inline int first() {
      return next[_head];
    }

    inline void add(const int elt)
    {
      next[_tail] = elt;
      _tail = elt;
      next[elt] = _head;
    }

    /// Pop the first event of the queue
    inline int pop()
    {
      int elt=next[_head];
      next[_head] = next[elt];
      if( next[_head] == _head ) _tail = _head;

      return elt;
    }

    /// Is the queue empty?
    inline bool empty() const
    {
      return _head == _tail;
    }

    /// Clear all elements of the queue
    inline void clear()
    {
      _tail = _head;
    }
    //@}

    /*!@name Miscellanous*/
    //@{ 
//     std::string getString() const
//     {
//       int i;
//       std::string return_str = "(";
//       i = next[_head];
//       if(i != _head) {
// 	return_str += (toString(i));
// 	i = next[i];
// 	while( i != _head ) {
// 	  return_str += (" "+toString(i));
// 	  i = next[i];
// 	}
//       }
//       return_str += ")";
//       return return_str;
//     }
    std::ostream& display(std::ostream& os) const {
      int i;
      os << "(";
      i = next[_head];
      if(i != _head) {
	os << i;
	i = next[i];
	while( i != _head ) {
	  os << " " << i;
	  i = next[i];
	}
      }
      os << ")";
      return os;
    }
    //@}

  };


//   template < class DATA_TYPE > 
//   std::string toString(const Vector< DATA_TYPE >& x) {
//     std::ostringstream os;
//     os << x;
//     return os.str();
//   }
  
//   template < class DATA_TYPE > 
//   std::string toString(const Stack< DATA_TYPE >& x) {
//     std::ostringstream os;
//     os << x;
//     return os.str();
//   }
  
//   template < class DATA_TYPE > 
//   std::string toString(const Array< DATA_TYPE >& x) {
//     return x.getString();
//   }

//   template < class DATA_TYPE > 
//   std::string toString(const Node< DATA_TYPE >& x) {
//     return x.getString();
//   }

//   template < class DATA_TYPE, int NUM_HEAD > 
//   std::string toString(const MultiList< DATA_TYPE, NUM_HEAD >& x) {
//     return x.getString();
//   }

//   template< class WORD_TYPE, class FLOAT_TYPE >
//   std::string toString(const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
//     return x.getString();
//   }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const IdStack< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE, int NUM_HEAD > 
  std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >& x) {
    return x.display(os);
  }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >& x) {
    return x.display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Vector< DATA_TYPE >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Stack< DATA_TYPE >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Array< DATA_TYPE >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE > 
  std::ostream& operator<< (std::ostream& os, const Node< DATA_TYPE >* x) {
    return x->display(os);
  }

  template < class DATA_TYPE, int NUM_HEAD > 
  std::ostream& operator<< (std::ostream& os, const MultiList< DATA_TYPE, NUM_HEAD >* x) {
    return x->display(os);
  }

  template< class WORD_TYPE, class FLOAT_TYPE >
  std::ostream& operator<< (std::ostream& os, const Bitset< WORD_TYPE, FLOAT_TYPE >* x) {
    return x->display(os);
  }

}


#endif // __STRUCTURE_HPP
