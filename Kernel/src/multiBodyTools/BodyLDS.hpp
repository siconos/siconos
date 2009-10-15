#ifndef BodyLDS_hpp
#define BodyLDS_hpp

class BodyLDS : public LagrangianDS,
  public boost::enable_shared_from_this<BodyLDS>
{
public:

  virtual selfHash(SP::SpaceFilter) = 0;

  virtual selfFindInteractions(SP::SpaceFilter) = 0;

  ACCEPT_VISITORS();
}

#endif
