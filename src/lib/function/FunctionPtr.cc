#include <config.h>
#include <function/FunctionPtr.h>

using std::string;

namespace jags {

FunctionPtr::FunctionPtr()
    : lfunc(nullptr), sfunc(nullptr), vfunc(nullptr), afunc(nullptr)
{}

FunctionPtr::FunctionPtr(ScalarFunction const *sf)
    : lfunc(nullptr), sfunc(sf), vfunc(nullptr), afunc(nullptr)
{}

FunctionPtr::FunctionPtr(VectorFunction const *vf)
    : lfunc(nullptr), sfunc(nullptr), vfunc(vf), afunc(nullptr)
{}

FunctionPtr::FunctionPtr(ArrayFunction const *af)
    : lfunc(nullptr), sfunc(nullptr), vfunc(nullptr), afunc(af)
{}

FunctionPtr::FunctionPtr(LinkFunction const *lf)
    : lfunc(lf), sfunc(nullptr), vfunc(nullptr), afunc(nullptr)
{}

string const &FunctionPtr::name() const
{
    static const string nullstring;
    if (lfunc)
	return lfunc->name();
    else if (sfunc)
	return sfunc->name();
    else if (vfunc)
	return vfunc->name();
    else if (afunc)
	return afunc->name();
    else
	return nullstring;
}

LinkFunction const * LINK(FunctionPtr const &p) { 
    return p.lfunc; 
}

ScalarFunction const * SCALAR(FunctionPtr const &p) { 
    return p.sfunc; 
}

VectorFunction const * VECTOR(FunctionPtr const &p) { 
    return p.vfunc; 
}

ArrayFunction const * ARRAY(FunctionPtr const &p) { 
    return p.afunc; 
}

Function const * FUNC(FunctionPtr const &p) 
{
    if (p.lfunc)
	return p.lfunc;
    else if (p.sfunc)
	return p.sfunc;
    else if (p.vfunc)
	return p.vfunc;
    else if (p.afunc)
	return p.afunc;
    else
	return nullptr;
}

bool FunctionPtr::operator==(FunctionPtr const &rhs) const
{
    return (lfunc == rhs.lfunc && sfunc==rhs.sfunc && vfunc==rhs.vfunc &&
	    afunc==rhs.afunc);
}
	

bool isNULL(FunctionPtr const &p)
{
    return (p.lfunc==nullptr && p.sfunc == nullptr && p.vfunc == nullptr && p.afunc == nullptr);
}

} //namespace jags
