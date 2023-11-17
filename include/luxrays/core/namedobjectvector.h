/***************************************************************************
 * Copyright 1998-2020 by authors (see AUTHORS.txt)                        *
 *                                                                         *
 *   This file is part of LuxCoreRender.                                   *
 *                                                                         *
 * Licensed under the Apache License, Version 2.0 (the "License");         *
 * you may not use this file except in compliance with the License.        *
 * You may obtain a copy of the License at                                 *
 *                                                                         *
 *     http://www.apache.org/licenses/LICENSE-2.0                          *
 *                                                                         *
 * Unless required by applicable law or agreed to in writing, software     *
 * distributed under the License is distributed on an "AS IS" BASIS,       *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 * See the License for the specific language governing permissions and     *
 * limitations under the License.                                          *
 ***************************************************************************/

#ifndef _LUXRAYS_NAMEDOBJECTVECTOR_H
#define	_LUXRAYS_NAMEDOBJECTVECTOR_H

#include <string>
#include <vector>
#include <ostream>

#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>

#include "luxrays/luxrays.h"
#include "luxrays/core/namedobject.h"

namespace luxrays {

	// no npreserva l'ordine originale in caso di delete
	template <class T>
	class FastVector
	{
	public:
		struct iterator : std::iterator<std::forward_iterator_tag, T>
		{
			public:
			iterator(T * p, int i) : pnt_(p), i_(i) {}
			iterator() {}

			T& operator*() {
				return  * (pnt_+i_);
			}

			iterator& operator++() {
				++i_;
				return *this;
			}

			bool operator==(const iterator& r) const {
				return i_ == r.i_;
			}
			bool operator!=(const iterator& r) const {
				return i_ != r.i_;
			}

			T* pnt_;
			int i_;
		};

		// Constructors/Destructor
		FastVector()
			:_size(0), _elements(new T[32]), _space(32)
		{}

		FastVector(int s)
			: _size(s), _elements(new T[s]), _space(s)
		{
			for (int index = 0; index < _size; ++index)
				_elements[index] = T();
		}

		~FastVector()
		{
			delete[] _elements;
		}

		inline int size() const noexcept {
			return _size;
		}

		T& operator[](int pos) {
			return _elements[pos];
		}

		const T& operator[](int pos) const {
			return _elements[pos];
		}

		//iterator* begin() noexcept;
		iterator begin() const noexcept {
			return FastVector<T>::iterator(_elements,0);
		}

		//iterator* end() noexcept;
		iterator end() const noexcept {
			return FastVector<T>::iterator(_elements, _size-1);
		}

		void push_back(const T& value)
		{
			if (_space == 0)
				reserve(8);
			else if (_size == _space)
				reserve(2 * _space);

			_elements[_size] = value;

			++_size;
		}
		void push_back(T&& value) {
			if (_space == 0)
				reserve(8);
			else if (_size == _space)
				reserve(2 * _space);

			_elements[_size] = value;

			++_size;
		}

		void eraseAt(int index) {
			if (index < _size - 1)
				_elements[index] = _elements[_size-1];
			_size--;
		}
		void pop( ) {
			_size--;
		}


		inline void reserve(int newalloc)
		{
			if (newalloc <= _space) return;

			T* p = new T[newalloc];

			for (int i = 0; i < _size; ++i)
				p[i] = _elements[i];

			delete[] _elements;

			_elements = p;

			_space = newalloc;
		}

	private:
		size_t	_size;		// Number of elements in Vector
		T* _elements;	// Pointer to first element of Vector
		size_t	_space;		// Total space used by Vector including
	};

// ==================================================================

class NamedObjectVector {
public:
	NamedObjectVector();
	virtual ~NamedObjectVector();

	NamedObject *DefineObj(NamedObject *newObj);
	bool IsObjDefined(const std::string &name) const;

	const NamedObject *GetObj(const std::string &name) const;
	NamedObject *GetObj(const std::string &name);
	const NamedObject *GetObj(const u_int index) const;
	NamedObject *GetObj(const u_int index);

	u_int GetIndex(const std::string &name) const;
	u_int GetIndex(const NamedObject *o) const;

	const std::string &GetName(const u_int index) const;
	const std::string &GetName(const NamedObject *o) const;

	u_int GetSize()const;
	void GetNames(std::vector<std::string> &names) const;
	std::vector<NamedObject *> &GetObjs();

	void DeleteObj(const std::string &name);
	void DeleteObjs(const std::vector<std::string> &names);
	
	std::string ToString() const;

private:
	typedef boost::bimap<boost::bimaps::unordered_set_of<std::string>,
			boost::bimaps::unordered_set_of<u_int> > Name2IndexType;
	typedef boost::bimap<boost::bimaps::unordered_set_of<u_int>,
			boost::bimaps::unordered_set_of<const NamedObject *> > Index2ObjType;

	//std::vector<NamedObject *> objs;
	std::vector<NamedObject*> objs;

	Name2IndexType name2index;
	Index2ObjType index2obj;
};

}

#endif	/* _LUXRAYS_NAMEDOBJECTVECTOR_H */
