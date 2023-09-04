#ifndef CL3DModel_includeguard
#define CL3DModel_includeguard


/*
Copyright 2023 Christian Linhart

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <memory>
#include <vector>
#include <unordered_set>
#include <list>
#include <eigen3/Eigen/Dense>

namespace CL3DModel {

class EdgeCycle;
class Edge;

class Vertex
{
public:
    Vertex()
    : _coords(0.0, 0.0, 0.0)
    {}

    Vertex(const Eigen::Vector3d& coords)
    : _coords(coords)
    {}

private:
    Eigen::Vector3d _coords;
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

    inline const QuarterEdgeRef& neighbor() const;
    QuarterEdgeRef otherEnd() const
    {
        return QuarterEdgeRef(_edge, _endPointIdx ^ 1, _edgeCycleIdx);
    }

    QuarterEdgeRef otherSide() const
    {
        return QuarterEdgeRef(_edge, _endPointIdx, _edgeCycleIdx ^ 1);
    }

    //getters
    inline Vertex* endPoint() const;
    inline Vertex* otherEndPoint() const;

    inline EdgeCycle* edgeCycle() const;
    inline EdgeCycle* otherSideEdgeCycle() const;

    Edge* edge() const { return _edge; }

    //operations
    

private:
    Edge* _edge;
    unsigned int _endPointIdx : 1;
    unsigned int _edgeCycleIdx : 1;
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


    Edge( Vertex* endPoint0, Vertex* endPoint1, const std::shared_ptr<Curve>& curve = std::shared_ptr<Curve>() )
    : _curve(curve)
    , _endPoints{ endPoint0, endPoint1 }
    , _edgeCycles{ nullptr, nullptr }
    {}
   
    operator bool() const { return _endPoints[0] && _endPoints[1]; }
    bool isNull() const { return !this->operator bool(); }



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

    EdgeCycle(const QuarterEdgeRef& oneEdge)
    : _oneEdge(oneEdge)
    {}

private:
    //one quarter-edge of the edge cycle. The others can be found through navigation to neighbor and otherEndpoint alternately
    QuarterEdgeRef _oneEdge;
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
    Plane( const Eigen::Vector3d& normalVector, double d );

    //distance from plane in units of length of _normalVector, sign indicates which side of the plane
    virtual double signedDistance( const Eigen::Vector3d& point )
    {
        return point.dot(_normalVector) + _d;
    }

private:
    Eigen::Vector3d _normalVector; //for signedDistance to yield the actual distance, _normalVector has to have length 1
    double _d;  // so that the plane equation is p.dot(_normalVector) + _d == 0
};

template <typename TEdgeCycleType = EdgeCycle>
class Facet
{
    typedef TEdgeCycleType EdgeCycleType;

    Facet( const std::shared_ptr<Surface>& surface )
    : _surface(surface)
    {}

    virtual~ Facet() {}

    EdgeCycleType *createEdgeCycle(const QuarterEdgeRef& oneEdge)
    {
        return &_edgeCycles.emplace_back(oneEdge);
    }

private:
    std::list<EdgeCycleType> _edgeCycles;
    std::shared_ptr<Surface> _surface;
};


template <typename TVertexType = Vertex, typename TEdgeType = Edge, typename TFacetType = Facet<> >
class Volume
{
public:

    typedef TVertexType VertexType;
    typedef TEdgeType EdgeType;
    typedef TFacetType FacetType;

    VertexType* createVertex(const Eigen::Vector3d& coords) const
    {
        return &_vertices.emplace_back(coords);
    }

    EdgeType* createEdge(Vertex* endPoint0, Vertex* endPoint1, const std::shared_ptr<Curve> &curve = std::shared_ptr<Curve>()) const
    {
        return &_edges.emplace_back(endPoint0, endPoint1, curve);
    }

    FacetType* createFacet(const std::shared_ptr<Surface>& surface) const
    {
        return &_edges.emplace_back(surface);
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
    using VolumeType::VertexType;
    using VolumeType::EdgeType;
    using VolumeType::FacetType;

private:
    std::unordered_set<std::shared_ptr<VolumeType>> _volumes;
};


//inline functions outside of class bodies
inline const QuarterEdgeRef& QuarterEdgeRef::neighbor() const
{
    return _edge->_neighbors[_endPointIdx][_edgeCycleIdx];
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

} //end namespace CL3DModel

#endif
