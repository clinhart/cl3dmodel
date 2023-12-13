#ifndef CL3DModel_includeguard
#define CL3DModel_includeguard


/*
Copyright 2023 Christian Linhart

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#define CL3DModel_MajorVersion 0
#define CL3DModel_MinorVersion 1
#define CL3DModel_PatchVersion 0


#include <memory>
#include <vector>
#include <set>
#include <list>
#include <Eigen/Dense>

namespace CL3DModel {

class EdgeCycle;
class Edge;
class Vertex;
class Facet;

template<typename TVertexType = Vertex, typename TEdgeType = Edge, typename TFacetType = Facet>
struct ModelTypes
{
	typedef TVertexType Vertex;
	typedef TEdgeType Edge;
	typedef TFacetType Facet;
};

template <typename TModelTypes>
class Volume;

class VolumeBase
{
public:
	virtual ~VolumeBase() {}
};

class QuarterEdgeRef
{
public:

	QuarterEdgeRef()
	: _edge(nullptr)
	, _vertexIdx(0)
	, _edgeCycleIdx(0)
	{}

	QuarterEdgeRef(Edge* edge, unsigned int vertexIdx, unsigned int edgeCycleIdx)
	: _edge(edge)
	, _vertexIdx(vertexIdx)
	, _edgeCycleIdx(edgeCycleIdx)
	{}

	bool isNull() const { return ( _edge == nullptr ); }

	operator bool() const { return ( _edge != nullptr ); }

	//navigation
	void moveToOtherEnd()
	{
		_vertexIdx ^= 1;
	}

	void moveToOtherSide()
	{
		_edgeCycleIdx ^= 1;
	}

	void moveToNeighbor()
	{
		*this = neighbor();
	}

	void moveToNeighborAtOtherEnd()
	{
		*this = neighborAtOtherEnd();
	}

	[[deprecated]] void moveToNeighborAtOtherEndPoint()
	{
		return moveToNeighborAtOtherEnd();
	}

	void moveToNeighborAtOtherSide()
	{
		*this = neighborAtOtherSide();
	}

	inline const QuarterEdgeRef& neighbor() const;
	inline const QuarterEdgeRef& neighborAtOtherEnd() const;
	[[deprecated]] inline const QuarterEdgeRef& neighborAtOtherEndPoint() const {
		return neighborAtOtherEnd();
	}
	inline const QuarterEdgeRef& neighborAtOtherSide() const;

	QuarterEdgeRef otherEnd() const
	{
		return QuarterEdgeRef(_edge, _vertexIdx ^ 1, _edgeCycleIdx);
	}

	QuarterEdgeRef otherSide() const
	{
		return QuarterEdgeRef(_edge, _vertexIdx, _edgeCycleIdx ^ 1);
	}

	//linking
	inline void connectToNeighbor(QuarterEdgeRef neighbor) const;
	inline void disconnectFromNeighbor() const;

	//getters
	inline Vertex* getVertex() const;
	inline Vertex* getOtherVertex() const;

	inline EdgeCycle* getEdgeCycle() const;
	inline EdgeCycle* getOtherSideEdgeCycle() const;

	Edge* getEdge() const { return _edge; }

	unsigned int getVertexIdx() const { return _vertexIdx; }
	unsigned int getEdgeCycleIdx() const { return _edgeCycleIdx; }

	//comparisons
	bool operator==( const QuarterEdgeRef& other ) const
	{
		return
			( this->_edge == other._edge )
			&&
			( this->_vertexIdx == other._vertexIdx )
			&&
			( this->_edgeCycleIdx == other._edgeCycleIdx )
		;
	}

	bool operator!=( const QuarterEdgeRef& other ) const
	{
		return !operator==(other);
	}

	//operations
	inline void setVertex(Vertex* vertex);

	inline QuarterEdgeRef mergeWithNeighbor(const std::function<void(Edge*)>& edgeDeleter = [](Edge*) {});

	//deprecated
	//use getVertex or getOtherVertex instead
	[[deprecated]] inline Vertex* endPoint() const { return getVertex(); }
	[[deprecated]] inline Vertex* otherEndPoint() const { return getOtherVertex(); };

	//use getEdgeCycle or getOtherSideEdgeCycle instead
	[[deprecated]] inline EdgeCycle* edgeCycle() const { return getEdgeCycle(); }
	[[deprecated]] inline EdgeCycle* otherSideEdgeCycle() const { return getOtherSideEdgeCycle();  }

	//use getEdge instead
	[[deprecated]] Edge* edge() const { return _edge; }

private:
	friend class EdgeCycle;

	inline QuarterEdgeRef& neighborWritableRef() const;
	inline QuarterEdgeRef* neighborWritablePtr() const;

	inline void setEdgeCycle( EdgeCycle* edgeCycle );

	Edge* _edge;
	unsigned int _vertexIdx : 1;
	unsigned int _edgeCycleIdx : 1;
};

class Vertex
{
public:
	Vertex()
		: _coords{ 0.0, 0.0, 0.0 }
	{}

	Vertex(const Eigen::Vector3d& coords)
	: _coords(coords)
	{}

	Vertex(Eigen::Vector3d&& coords)
		: _coords(std::move(coords))
	{}

	virtual ~Vertex() {}

	const Eigen::Vector3d& getCoords() const { return _coords; }
	const QuarterEdgeRef& getOneEdge() const { return _oneEdge; }

	void indicateNewEdge( const QuarterEdgeRef& edge )
	{
		if (_oneEdge.isNull()) {
			_oneEdge = edge;
		}
	}

	void indicateBeforeRemoveEdge(const QuarterEdgeRef& edge)
	{
		if (_oneEdge == edge || _oneEdge == edge.otherSide()) {
			_oneEdge.moveToNeighbor();
			if (_oneEdge) {
				_oneEdge.moveToOtherSide();

				if (_oneEdge == edge || _oneEdge == edge.otherSide()) {
					_oneEdge = QuarterEdgeRef();
				}
			}
		}
	}

private:
    //quarter-edge reference to one of the quarter-edges that is at this vertex.
    //The others can be found through navigation to neighbor and otherEndpoint alternately
    QuarterEdgeRef _oneEdge;

    Eigen::Vector3d _coords;
};

//all curve objects are immutable, i.e. created by the constructor and then not changed
//this is needed so that multiple edges can share the same curve
class Curve
{
public:
    virtual ~Curve() {}
};

#ifdef CL3DMODEL_CHECK_SAME_COORDS
inline void SameCoords()
{

}
#endif
class Edge
{
public:
	Edge()
		: _curve()
		, _vertices{ nullptr, nullptr }
		, _edgeCycles{ nullptr, nullptr }
	{}

	virtual ~Edge() {}

	Edge(Vertex* vertex0, Vertex* vertex1, const std::shared_ptr<Curve>& curve = std::shared_ptr<Curve>())
		: _curve(curve)
		, _vertices{ vertex0, vertex1 }
		, _edgeCycles{ nullptr, nullptr }
	{
		for (int vertexIdx = 0; vertexIdx < 2; vertexIdx++) {
			if (_vertices[vertexIdx]) {
				_vertices[vertexIdx]->indicateNewEdge(
					QuarterEdgeRef(this, vertexIdx, 0)
				);
			}
		}
	}

	operator bool() const { return ( _vertices[0] || _vertices[1] || _curve ); }
	bool isNull() const { return !this->operator bool(); }

	//idx == 0: start-point, idx == 1: end-point
	Vertex* getVertex(int idx) const { assert((idx == 0) || (idx == 1)); return _vertices[idx]; }

	EdgeCycle* getEdgeCycle(int idx) const { assert((idx == 0) || (idx == 1)); return _edgeCycles[idx]; }

	std::shared_ptr<Curve> getCurve() const { return _curve; }

	//change edge 
	//set a vertex. This is only allowed if the edge has no neighors at that vertex
	void setVertex(int vertexIdx, Vertex* vertex) {
		assert(_neighbors[vertexIdx][0].isNull());
		assert(_neighbors[vertexIdx][1].isNull());
		assert((vertexIdx == 0) || (vertexIdx == 1));
		_vertices[vertexIdx] = vertex;
#ifdef CL3DMODEL_CHECK_SAME_COORDS
		checkSameCoords();
#endif
	}

	//deprecated
	//use getVertex instead
	[[deprecated]] Vertex* getEndPoint(int idx) const { return getVertex(idx); }

#ifdef CL3DMODEL_CHECK_SAME_COORDS
	void checkSameCoords() {
		if (_vertices[0] && _vertices[1] && (_vertices[0]->getCoords() - _vertices[1]->getCoords()).squaredNorm() < 0.0001) {
			SameCoords();
		}
	}
#endif

private:
    friend class QuarterEdgeRef;
	template <typename TModelTypes> friend class Volume;

	//function
	// 
	//set a vertex ( use with caution! this operation alone will render the model inconsistent! That's why it's a private function)
	void setVertexUnchecked(int vertexIdx, Vertex* vertex) {
		assert((vertexIdx == 0) || (vertexIdx == 1));
		_vertices[vertexIdx] = vertex;
	}

	//data
    std::shared_ptr<Curve> _curve; //if _curve is null, then the edge is a straight line
    Vertex* _vertices[2];
    EdgeCycle* _edgeCycles[2];

    QuarterEdgeRef _neighbors[2][2]; //neighbors[vertexIdx][edgeCycleIdx]
};

class EdgeCycle
{
public:

	EdgeCycle(const QuarterEdgeRef& oneEdge, Facet* facet)
	: _oneEdge(oneEdge)
	, _facet(facet)
	{
		//set the pointers from edges to this edgeCycle
		if (oneEdge) {
			QuarterEdgeRef qedge = oneEdge;
			do {
				qedge.setEdgeCycle(this);
				qedge.moveToNeighborAtOtherEnd();
			} while( qedge && ( qedge != oneEdge ) );
		}
	}

	virtual ~EdgeCycle() {}

	const QuarterEdgeRef& getOneEdge() const { return _oneEdge; }
	Facet * getFacet() const { return _facet; }

	void flipOrientation()
	{
		_oneEdge.moveToOtherEnd();
	}

private: friend class QuarterEdgeRef;

	void setOneEdge(const QuarterEdgeRef& oneEdge) { _oneEdge = oneEdge; }
private:
	//one quarter-edge of the edge cycle. The others can be found through navigation to neighbor and otherEndpoint alternately
	QuarterEdgeRef _oneEdge;
	Facet *_facet;
};

//all surface objects are immutable, i.e. created by the constructor and then not changed
//this is needed so that multiple facets can share the same surface
class Surface
{
public:
    virtual ~Surface() {}

    //signed distance from surface. Sign is relevant for a true volume model: positive number is outside of volume, negative number is inside of volume
    virtual double signedDistance( const Eigen::Vector3d& point ) = 0;
};

class Plane : public Surface
{

public:
    //so that the plane equation is p.dot(normalVector) + d == 0
	Plane(const Eigen::Vector3d& normalVector, double d)
		: _normalVector(normalVector)
		, _d(d)
	{
	}

    //distance from plane in units of length of _normalVector, sign indicates which side of the plane
    virtual double signedDistance( const Eigen::Vector3d& point )
    {
        return point.dot(_normalVector) + _d;
    }

	const Eigen::Vector3d& getNormalVector() const { return _normalVector; }
	double getD() const { return _d; }

private:
    Eigen::Vector3d _normalVector; //for signedDistance to yield the actual distance, _normalVector has to have length 1
    double _d;  // so that the plane equation is p.dot(_normalVector) + _d == 0
};

class Facet
{
public:
	typedef typename std::list<EdgeCycle>::iterator EdgeCycleIterator;
	typedef typename std::list<EdgeCycle>::const_iterator ConstEdgeCycleIterator;

    Facet( const std::shared_ptr<Surface>& surface )
    : _surface(surface)
    {}

    virtual ~Facet() {}

    EdgeCycle *createEdgeCycle(const QuarterEdgeRef& oneEdge)
    {
        return &_edgeCycles.emplace_back(oneEdge, this);
    }

	EdgeCycleIterator begin() { return _edgeCycles.begin(); }
	EdgeCycleIterator end() { return _edgeCycles.end(); }

	ConstEdgeCycleIterator cbegin() const { return _edgeCycles.cbegin(); }
	ConstEdgeCycleIterator cend() const { return _edgeCycles.cend(); }

	ConstEdgeCycleIterator begin() const { return _edgeCycles.cbegin(); }
	ConstEdgeCycleIterator end() const { return _edgeCycles.cend(); }

	size_t getEdgeCycleCount() const { return _edgeCycles.size(); }

	const std::shared_ptr<Surface>& getSurface() const { return _surface; }

	template<typename VolumeType> VolumeType* getVolume() const
	{
		return dynamic_cast<VolumeType*>(_volume);
	}

private:
	template <typename TModelTypes> friend class Volume;

	void setVolume(VolumeBase* volume) { _volume = volume; }

private:
    std::list<EdgeCycle> _edgeCycles;
    std::shared_ptr<Surface> _surface;
	VolumeBase* _volume;
};


template <typename TModelTypes = ModelTypes<> >
class Volume : public VolumeBase
{
public:

	typedef typename TModelTypes::Vertex VertexType;
	typedef typename TModelTypes::Edge EdgeType;
	typedef typename TModelTypes::Facet FacetType;

	class StorageEdgeType : public EdgeType {
		using EdgeType::EdgeType;

		friend class Volume;

	private:
		typename std::list<StorageEdgeType>::iterator _iterator;
	};

	typedef typename std::list<VertexType>::iterator VertexIterator;
	typedef typename std::list<StorageEdgeType>::iterator EdgeIterator;
	typedef typename std::list<FacetType>::iterator FacetIterator;

	typedef typename std::list<VertexType>::const_iterator ConstVertexIterator;
	typedef typename std::list<StorageEdgeType>::const_iterator ConstEdgeIterator;
	typedef typename std::list<FacetType>::const_iterator ConstFacetIterator;

	template<typename... ArgTypes>
	VertexType* createVertex(ArgTypes... args)
	{
		return &_vertices.emplace_back(args...);
	}

	template<typename... ArgTypes>
	EdgeType* createEdge(ArgTypes... args)
	{
		StorageEdgeType* edge = &_edges.emplace_back(args...);
		edge->_iterator = std::prev(_edges.end());
		return edge;
	}

	void deleteEdge(EdgeType* edge)
	{
		if (edge) {
			StorageEdgeType* storageEdge = static_cast<StorageEdgeType*>(edge);
			auto itEdge = storageEdge->_iterator;
			if (itEdge != _edges.end()) {
				storageEdge->_iterator = _edges.end();
				_edges.erase(itEdge);
			}
		}
	}

	template<typename... ArgTypes>
	FacetType* createFacet(ArgTypes... args)
	{
		FacetType* facet = &_facets.emplace_back(args...);
		facet->setVolume(this);
		return facet;
	}

	//iterate over content
	ConstVertexIterator cbeginVertices() const { return _vertices.cbegin(); }
	ConstVertexIterator cendVertices() const { return _vertices.cend(); }

	ConstEdgeIterator cbeginEdges() const { return _edges.cbegin(); }
	ConstEdgeIterator cendEdges() const { return _edges.cend(); }

	ConstFacetIterator cbeginFacets() const { return _facets.cbegin(); }
	ConstFacetIterator cendFacets() const { return _facets.cend(); }

	VertexIterator beginVertices() { return _vertices.begin(); }
	VertexIterator endVertices() { return _vertices.end(); }

	EdgeIterator beginEdges() { return _edges.begin(); }
	EdgeIterator endEdges() { return _edges.end(); }

	FacetIterator beginFacets() { return _facets.begin(); }
	FacetIterator endFacets() { return _facets.end(); }

	ConstVertexIterator beginVertices() const { return _vertices.cbegin(); }
	ConstVertexIterator endVertices() const { return _vertices.cend(); }

	ConstEdgeIterator beginEdges() const { return _edges.cbegin(); }
	ConstEdgeIterator endEdges() const { return _edges.cend(); }

	ConstFacetIterator beginFacets() const { return _facets.cbegin(); }
	ConstFacetIterator endFacets() const { return _facets.cend(); }

	//get quantities
	size_t getVertexCount() const { return _vertices.size(); }
	size_t getEdgeCount() const { return _edges.size(); }
	size_t getFacetCount() const { return _facets.size(); }

	//operations
	//split edge at the given vertex in two edges
	//the edge should be of this volume
	//returns a pair of the two resuling edges: the first part has the same vertex at idx=0
	//one of the parts may have the same pointer value as the given edge, but this is not guaranteed
	template<typename... AdditionalCreateEdgeArgsTypes>
	std::pair<EdgeType*, EdgeType*> splitEdge(EdgeType* edge, VertexType* vertex, AdditionalCreateEdgeArgsTypes... additionalCreateEdgeArgs)
	{
		EdgeType* newEdge =
			this->createEdge(
				vertex,
				edge->getVertex(1),
				edge->getCurve(),
				additionalCreateEdgeArgs...
			);

		QuarterEdgeRef oldEndsOfEdge[2];
		QuarterEdgeRef oldNeighborsAtEndOfEdge[2];

		//get neighbors at end of edge
		for (int cycleIdx = 0; cycleIdx < 2; cycleIdx++) {
			oldEndsOfEdge[cycleIdx] = QuarterEdgeRef(edge, 1, cycleIdx);
			oldNeighborsAtEndOfEdge[cycleIdx] = oldEndsOfEdge[cycleIdx].neighbor();
		}

		//connect new edge with those neighbors
		for (int cycleIdx = 0; cycleIdx < 2; cycleIdx++) {
			if (oldNeighborsAtEndOfEdge[cycleIdx]) {
				QuarterEdgeRef(newEdge, 1, cycleIdx).connectToNeighbor(
					oldNeighborsAtEndOfEdge[cycleIdx]
				);
			}
		}

		//change end vertex of edge
		edge->setVertexUnchecked(1, vertex);

		//connect edge and new edge as neighbors for both cycles
		for (int cycleIdx = 0; cycleIdx < 2; cycleIdx++) {
			QuarterEdgeRef(edge, 1, cycleIdx).connectToNeighbor(
				QuarterEdgeRef(newEdge, 0, cycleIdx)
			);
		}

#ifdef CL3DMODEL_CHECK_SAME_COORDS
		edge->checkSameCoords();
		newEdge->checkSameCoords();
#endif
		//return the two parts
		return std::pair<EdgeType*, EdgeType*>(edge, newEdge);
	}

	//same as above but work with QuarterEdgeRefs
	//the first of the resulting pair has the same vertex as the gioven qEdge.
	template<typename... AdditionalCreateEdgeArgsTypes>
	std::pair<QuarterEdgeRef, QuarterEdgeRef> splitQEdge(const QuarterEdgeRef& qEdge, VertexType* vertex, AdditionalCreateEdgeArgsTypes... additionalCreateEdgeArgs)
	{
		std::pair<EdgeType*, EdgeType*> splitEdgeResult =
			splitEdge(static_cast<EdgeType*>(qEdge.getEdge()), vertex, additionalCreateEdgeArgs...);
		if (qEdge.getVertexIdx() == 0) {
			return std::pair<QuarterEdgeRef, QuarterEdgeRef>(
				QuarterEdgeRef(splitEdgeResult.first, 0, qEdge.getEdgeCycleIdx()),
				QuarterEdgeRef(splitEdgeResult.second, 0, qEdge.getEdgeCycleIdx())
			);
		} else {
			return std::pair<QuarterEdgeRef, QuarterEdgeRef>(
				QuarterEdgeRef(splitEdgeResult.second, 1, qEdge.getEdgeCycleIdx()),
				QuarterEdgeRef(splitEdgeResult.first, 1, qEdge.getEdgeCycleIdx())
			);
		}
	}


	QuarterEdgeRef mergeWithNeighbor(QuarterEdgeRef qEdge)
	{
		return qEdge.mergeWithNeighbor([this](Edge* edge) { this->deleteEdge(static_cast<EdgeType*>(edge)); });
	}
private:
	template <typename TVolumeType>
	friend class Scene;

	typedef uint64_t InternalIdType;

	InternalIdType _internalId;


	//TODO use more efficient storage than std::list. Elements need to stay at the same storage address, so std::vector is not possible
	std::list<VertexType> _vertices;
	std::list<StorageEdgeType> _edges;
	std::list<FacetType> _facets;
};


template <typename TVolumeType = Volume<>>
class Scene
{
public:
	typedef TVolumeType VolumeType;
	typedef typename VolumeType::VertexType VertexType;
	typedef typename VolumeType::EdgeType EdgeType;
	typedef typename VolumeType::FacetType FacetType;

	typedef std::shared_ptr<VolumeType> VolumePointerType;

	typedef typename std::set<VolumePointerType>::const_iterator const_iterator;

	template<typename... ArgTypes>
	VolumePointerType createVolume(ArgTypes... args)
	{
		VolumePointerType volume = std::make_shared<VolumeType>(args...);
		volume->_internalId = _nextInternalVolumeId;
		++_nextInternalVolumeId;
		_volumes.insert(volume);
		return volume;
	}

	//iterate over volumes
	const_iterator begin() const { return _volumes.begin(); }
	const_iterator end() const { return _volumes.end(); }

	//get quantities
	size_t getVolumesCount() const { return _volumes.size(); }

private:
	typename TVolumeType::InternalIdType _nextInternalVolumeId = 1;
	typedef std::set<
		VolumePointerType,
		bool (*)(
			const VolumePointerType& vol1,
			const VolumePointerType& vol2
			)
	> VolumeSetType;

	static bool cmpVolumesByInternalId(
		const VolumePointerType& vol1,
		const VolumePointerType& vol2
	) {
		return vol1->_internalId < vol2->_internalId;
	}

	VolumeSetType _volumes{ cmpVolumesByInternalId };
};


//inline functions outside of class bodies
inline const QuarterEdgeRef& QuarterEdgeRef::neighbor() const
{
    return _edge->_neighbors[_vertexIdx][_edgeCycleIdx];
}

inline const QuarterEdgeRef& QuarterEdgeRef::neighborAtOtherEnd() const
{
    return _edge->_neighbors[_vertexIdx ^ 1][_edgeCycleIdx];
}

inline const QuarterEdgeRef& QuarterEdgeRef::neighborAtOtherSide() const
{
    return _edge->_neighbors[_vertexIdx][_edgeCycleIdx ^ 1];
}

inline QuarterEdgeRef& QuarterEdgeRef::neighborWritableRef() const
{
    return _edge->_neighbors[_vertexIdx][_edgeCycleIdx];
}

inline QuarterEdgeRef* QuarterEdgeRef::neighborWritablePtr() const
{
    return &(_edge->_neighbors[_vertexIdx][_edgeCycleIdx]);
}

inline void QuarterEdgeRef::connectToNeighbor(QuarterEdgeRef neighbor) const
{
    assert( !this->isNull() );
    assert( !neighbor.isNull() );
//    assert(this->getEdgeCycle() == nullptr || neighbor.getEdgeCycle() == nullptr || this->getEdgeCycle() == neighbor.getEdgeCycle() );
	assert( this->getVertex() == neighbor.getVertex() );
    QuarterEdgeRef* backLink = neighbor.neighborWritablePtr();
    if (! backLink->isNull() ) {
        backLink->disconnectFromNeighbor();
    }
    assert( backLink->isNull() );
    QuarterEdgeRef* forwardLink = this->neighborWritablePtr();
	if (!forwardLink->isNull()) {
		forwardLink->disconnectFromNeighbor();
	}

    *forwardLink = neighbor;
    *backLink = *this;

	EdgeCycle* thisEdgeCycle = this->getEdgeCycle();
	EdgeCycle* neighborEdgeCycle = neighbor.getEdgeCycle();
	if (thisEdgeCycle == nullptr) {
		if (neighborEdgeCycle) {
			this->_edge->_edgeCycles[_edgeCycleIdx] = neighborEdgeCycle;
		}
	} else if (neighborEdgeCycle == nullptr) {
		if (thisEdgeCycle) {
			neighbor._edge->_edgeCycles[_edgeCycleIdx] = thisEdgeCycle;
		}
	}
}

inline void QuarterEdgeRef::disconnectFromNeighbor() const
{
    assert( !this->isNull() );
    QuarterEdgeRef* forwardLink = this->neighborWritablePtr();
    if (forwardLink->_edge) {
        QuarterEdgeRef* backLink = forwardLink->neighborWritablePtr();

        assert(*backLink == *this);
        assert(forwardLink->_edge == backLink->neighbor()._edge);
        assert(backLink->_edge == forwardLink->neighbor()._edge);

        forwardLink->_edge = nullptr;
        backLink->_edge = nullptr;
    }
}

inline Vertex* QuarterEdgeRef::getVertex() const
{
    return _edge->_vertices[_vertexIdx];
}

inline Vertex* QuarterEdgeRef::getOtherVertex() const
{
    return _edge->_vertices[_vertexIdx ^ 1];
}

inline EdgeCycle* QuarterEdgeRef::getEdgeCycle() const
{
    return _edge->_edgeCycles[_edgeCycleIdx];
}

inline EdgeCycle* QuarterEdgeRef::getOtherSideEdgeCycle() const
{
    return _edge->_edgeCycles[_edgeCycleIdx ^ 1];
}

inline void QuarterEdgeRef::setEdgeCycle( EdgeCycle* edgeCycle )
{
	_edge->_edgeCycles[_edgeCycleIdx] = edgeCycle;
}

inline void QuarterEdgeRef::setVertex(Vertex* vertex)
{
	_edge->setVertex(_vertexIdx, vertex);
}


inline QuarterEdgeRef QuarterEdgeRef::mergeWithNeighbor(const std::function<void(Edge*)>& edgeDeleter)
{
	assert(neighbor());
	assert(neighbor().otherSide() == otherSide().neighbor());

	QuarterEdgeRef myNeighbor = neighbor();
	QuarterEdgeRef otherEndOfMyNeighbor = myNeighbor.otherEnd();
	Vertex* vertexOnOtherEndOfMyNeighbor = otherEndOfMyNeighbor.getVertex();
	QuarterEdgeRef nextNeighborOfMyNeighbor = otherEndOfMyNeighbor.neighbor();

	QuarterEdgeRef myOtherSide = otherSide();
	QuarterEdgeRef myNeighborOnOtherSide = myOtherSide.neighbor();
	QuarterEdgeRef otherEndOfMyNeighborOnOtherSide = myNeighborOnOtherSide.otherEnd();
	QuarterEdgeRef nextNeighborOfMyNeighborOnOtherSide = otherEndOfMyNeighborOnOtherSide.neighbor();

	//indicate to the vertex that myNeighbor will be deleted -> this removes a potential pointer from the vertex to myNeighbor
	Vertex* vertex = getVertex();
	vertex->indicateBeforeRemoveEdge(myNeighbor);

	//edgecycle stuff
	EdgeCycle* edgeCycleOfMyNeighbor = myNeighbor.getEdgeCycle();
	if (edgeCycleOfMyNeighbor && edgeCycleOfMyNeighbor->getOneEdge() == myNeighbor) {
		edgeCycleOfMyNeighbor->setOneEdge(*this);
	}

	EdgeCycle* edgeCycleOfMyNeighborOnOtherSide = myNeighborOnOtherSide.getEdgeCycle();
	if (edgeCycleOfMyNeighborOnOtherSide && edgeCycleOfMyNeighborOnOtherSide->getOneEdge() == myNeighborOnOtherSide) {
		edgeCycleOfMyNeighbor->setOneEdge(myOtherSide);
	}

	myNeighbor.setEdgeCycle(nullptr);
	myNeighborOnOtherSide.setEdgeCycle(nullptr);

	//disconnect edges
	this->disconnectFromNeighbor();
	myOtherSide.disconnectFromNeighbor();
	otherEndOfMyNeighbor.disconnectFromNeighbor();
	otherEndOfMyNeighborOnOtherSide.disconnectFromNeighbor();

	//delete the other edge
	edgeDeleter(myNeighbor.edge());

	//set vertex
	setVertex(vertexOnOtherEndOfMyNeighbor);

	//connect edges
	this->connectToNeighbor(nextNeighborOfMyNeighbor);
	myOtherSide.connectToNeighbor(nextNeighborOfMyNeighborOnOtherSide);

	return *this;
}

} //end namespace CL3DModel

#endif
