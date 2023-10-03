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
#include <unordered_set>
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


class QuarterEdgeRef
{
public:

	QuarterEdgeRef()
	: _edge(nullptr)
	, _endPointIdx(0)
	, _edgeCycleIdx(0)
	{}

	QuarterEdgeRef(Edge* edge, unsigned int endPointIdx, unsigned int edgeCycleIdx)
	: _edge(edge)
	, _endPointIdx(endPointIdx)
	, _edgeCycleIdx(edgeCycleIdx)
	{}

	bool isNull() const { return ( _edge == nullptr ); }

	operator bool() const { return ( _edge != nullptr ); }

	//navigation
	void moveToOtherEnd()
	{
		_endPointIdx ^= 1;
	}

	void moveToOtherSide()
	{
		_edgeCycleIdx ^= 1;
	}

	void moveToNeighbor()
	{
		*this = neighbor();
	}

	void moveToNeighborAtOtherEndPoint()
	{
		*this = neighborAtOtherEndPoint();
	}

	void moveToNeighborAtOtherSide()
	{
		*this = neighborAtOtherSide();
	}

	inline const QuarterEdgeRef& neighbor() const;
	inline const QuarterEdgeRef& neighborAtOtherEndPoint() const;
	inline const QuarterEdgeRef& neighborAtOtherSide() const;


	QuarterEdgeRef otherEnd() const
	{
		return QuarterEdgeRef(_edge, _endPointIdx ^ 1, _edgeCycleIdx);
	}

	QuarterEdgeRef otherSide() const
	{
		return QuarterEdgeRef(_edge, _endPointIdx, _edgeCycleIdx ^ 1);
	}

	//linking
	inline void connectToNeighbor(const QuarterEdgeRef& neighbor) const;
	inline void disconnectFromNeighbor() const;

	//getters
	inline Vertex* endPoint() const;
	inline Vertex* otherEndPoint() const;

	inline EdgeCycle* edgeCycle() const;
	inline EdgeCycle* otherSideEdgeCycle() const;

	Edge* edge() const { return _edge; }

	//comparisons
	bool operator==( const QuarterEdgeRef& other ) const
	{
		return
			( this->_edge == other._edge )
			&&
			( this->_endPointIdx == other._endPointIdx )
			&&
			( this->_edgeCycleIdx == other._edgeCycleIdx )
		;
	}

	bool operator!=( const QuarterEdgeRef& other ) const
	{
		return !operator==(other);
	}

	//operations


private:
	friend class EdgeCycle;

	inline QuarterEdgeRef& neighborWritableRef() const;
	inline QuarterEdgeRef* neighborWritablePtr() const;

	inline void setEdgeCycle( EdgeCycle* edgeCycle );

	Edge* _edge;
	unsigned int _endPointIdx : 1;
	unsigned int _edgeCycleIdx : 1;
};

class Vertex
{
public:
	Vertex()
	: _coords(0.0, 0.0, 0.0)
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

class Edge
{
public:
	Edge()
	: _curve()
	, _endPoints{ nullptr, nullptr }
	, _edgeCycles{ nullptr, nullptr }
	{}

	virtual ~Edge() {}

	Edge( Vertex* endPoint0, Vertex* endPoint1, const std::shared_ptr<Curve>& curve = std::shared_ptr<Curve>() )
	: _curve(curve)
	, _endPoints{ endPoint0, endPoint1 }
	, _edgeCycles{ nullptr, nullptr }
	{
		for( int endPointIdx = 0; endPointIdx < 2; endPointIdx++ ) {
			_endPoints[endPointIdx]->indicateNewEdge(
				QuarterEdgeRef(this, endPointIdx, 0)
			);
		}
	}

	operator bool() const { return _endPoints[0] && _endPoints[1]; }
	bool isNull() const { return !this->operator bool(); }

	//idx == 0: start-point, idx == 1: end-point
	Vertex* getEndPoint(int idx) const { return _endPoints[idx]; }

	EdgeCycle* getEdgeCycle(int idx) const { return _edgeCycles[idx]; }

private:
    friend class QuarterEdgeRef;

    std::shared_ptr<Curve> _curve; //if _curve is null, then the edge is a straight line
    Vertex* _endPoints[2];
    EdgeCycle* _edgeCycles[2];

    QuarterEdgeRef _neighbors[2][2]; //neighbors[endPointIdx][edgeCycleIdx]
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
				qedge.moveToNeighborAtOtherEndPoint();
			} while( qedge && ( qedge != oneEdge ) );
		}
	}

	virtual ~EdgeCycle() {}

	const QuarterEdgeRef& getOneEdge() const { return _oneEdge; }
	Facet * getFacet() const { return _facet; }

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

private:
    Eigen::Vector3d _normalVector; //for signedDistance to yield the actual distance, _normalVector has to have length 1
    double _d;  // so that the plane equation is p.dot(_normalVector) + _d == 0
};

class Facet
{
public:
    Facet( const std::shared_ptr<Surface>& surface )
    : _surface(surface)
    {}

    virtual ~Facet() {}

    EdgeCycle *createEdgeCycle(const QuarterEdgeRef& oneEdge)
    {
        return &_edgeCycles.emplace_back(oneEdge, this);
    }

private:
    std::list<EdgeCycle> _edgeCycles;
    std::shared_ptr<Surface> _surface;
};


template <typename TModelTypes = ModelTypes<> >
class Volume
{
public:

	typedef typename TModelTypes::Vertex VertexType;
	typedef typename TModelTypes::Edge EdgeType;
	typedef typename TModelTypes::Facet FacetType;

	template<typename... ArgTypes>
	VertexType* createVertex(ArgTypes... args)
	{
		return &_vertices.emplace_back(args...);
	}

	template<typename... ArgTypes>
	EdgeType* createEdge(ArgTypes... args)
	{
		return &_edges.emplace_back(args...);
	}

	template<typename... ArgTypes>
	FacetType* createFacet(ArgTypes... args)
	{
		return &_facets.emplace_back(args...);
	}

	//operations
	//split edge at the given vertex in two edges
	//the edge should be of this volume
	//TODO// void SplitEdge(EdgeType* edge, VertexType* vertex);
private:

	//TODO use more efficient storage than std::list. Elements need to stay at the same storage address, so std::vector is not possible
	std::list<VertexType> _vertices;
	std::list<EdgeType> _edges;
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

	template<typename... ArgTypes>
	VolumePointerType createVolume( ArgTypes... args )
	{
		return std::make_shared<VolumeType>(args...);
	}

private:
    std::unordered_set<std::shared_ptr<VolumeType>> _volumes;
};


//inline functions outside of class bodies
inline const QuarterEdgeRef& QuarterEdgeRef::neighbor() const
{
    return _edge->_neighbors[_endPointIdx][_edgeCycleIdx];
}

inline const QuarterEdgeRef& QuarterEdgeRef::neighborAtOtherEndPoint() const
{
    return _edge->_neighbors[_endPointIdx ^ 1][_edgeCycleIdx];
}

inline const QuarterEdgeRef& QuarterEdgeRef::neighborAtOtherSide() const
{
    return _edge->_neighbors[_endPointIdx][_edgeCycleIdx ^ 1];
}

inline QuarterEdgeRef& QuarterEdgeRef::neighborWritableRef() const
{
    return _edge->_neighbors[_endPointIdx][_edgeCycleIdx];
}

inline QuarterEdgeRef* QuarterEdgeRef::neighborWritablePtr() const
{
    return &(_edge->_neighbors[_endPointIdx][_edgeCycleIdx]);
}

inline void QuarterEdgeRef::connectToNeighbor(const QuarterEdgeRef& neighbor) const
{
    assert( !this->isNull() );
    assert( !neighbor.isNull() );
    assert( this->edgeCycle() == neighbor.edgeCycle() );
    QuarterEdgeRef* backLink = neighbor.neighborWritablePtr();
    if (! backLink->isNull() ) {
        backLink->disconnectFromNeighbor();
    }
    assert( backLink->isNull() );
    QuarterEdgeRef* forwardLink = this->neighborWritablePtr();

    *forwardLink = neighbor;
    *backLink = *this;
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

inline Vertex* QuarterEdgeRef::endPoint() const
{
    return _edge->_endPoints[_endPointIdx];
}

inline Vertex* QuarterEdgeRef::otherEndPoint() const
{
    return _edge->_endPoints[_endPointIdx ^ 1];
}

inline EdgeCycle* QuarterEdgeRef::edgeCycle() const
{
    return _edge->_edgeCycles[_edgeCycleIdx];
}

inline EdgeCycle* QuarterEdgeRef::otherSideEdgeCycle() const
{
    return _edge->_edgeCycles[_edgeCycleIdx ^ 1];
}

inline void QuarterEdgeRef::setEdgeCycle( EdgeCycle* edgeCycle )
{
	_edge->_edgeCycles[_edgeCycleIdx] = edgeCycle;
}


} //end namespace CL3DModel

#endif
