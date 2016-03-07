#undef PY_REGISTER
#undef PY_REGISTER_UNSHARED
#undef PY_REGISTER_WITHOUT_DIRECTOR

%define PY_REGISTER_UNSHARED(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR_REF_ONLY(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER_SIMPLEMATRIX(TYPE)
%include TYPE.hpp
%enddef

