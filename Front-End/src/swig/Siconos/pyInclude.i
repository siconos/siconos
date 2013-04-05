#undef PY_REGISTER
#undef PY_REGISTER_WITHOUT_DIRECTOR

%define PY_REGISTER(TYPE)
%include TYPE.hpp
%enddef

%define PY_REGISTER_WITHOUT_DIRECTOR(TYPE)
%include TYPE.hpp
%enddef
