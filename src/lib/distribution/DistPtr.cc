#include <config.h>
#include <distribution/DistPtr.h>

using std::string;

namespace jags {

DistPtr::DistPtr()
    : sdist(nullptr), vdist(nullptr), adist(nullptr)
{}

DistPtr::DistPtr(ScalarDist const *sf)
    : sdist(sf), vdist(nullptr), adist(nullptr)
{}

DistPtr::DistPtr(VectorDist const *vf)
    : sdist(nullptr), vdist(vf), adist(nullptr)
{}

DistPtr::DistPtr(ArrayDist const *af)
    : sdist(nullptr), vdist(nullptr), adist(af)
{}

string const &DistPtr::name() const
{
    static const string nullstring;
    if (sdist)
	return sdist->name();
    else if (vdist)
	return vdist->name();
    else if (adist)
	return adist->name();
    else
	return nullstring;
}

ScalarDist const * SCALAR(DistPtr const &p) { 
    return p.sdist; 
}

VectorDist const * VECTOR(DistPtr const &p) { 
    return p.vdist; 
}

ArrayDist const * ARRAY(DistPtr const &p) { 
    return p.adist; 
}

bool DistPtr::operator==(DistPtr const &rhs) const
{
    return (sdist==rhs.sdist && vdist==rhs.vdist && adist==rhs.adist);
}

bool isNULL(DistPtr const &p)
{
    return (p.sdist == nullptr && p.vdist == nullptr && p.adist == nullptr);
}

} //namespace jags
