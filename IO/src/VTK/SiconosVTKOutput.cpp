#include "SiconosVTKOutput.hpp"
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>

#include <SphereNEDS.hpp>

struct SiconosVTKOutput::_DataSetMaker : public SiconosVisitor
{
  SiconosVTKOutput& parent;

  _DataSetMaker(SiconosVTKOutput& p) : parent(p) {};

  void visit(const SphereNEDS& sphere)
  {
    vtkSmartPointer<vtkSphereSource> vtksphere =
      vtkSmartPointer<vtkSphereSource>::New();

    parent._sources.push_back(vtksphere);

    vtkSmartPointer<vtkPolyDataMapper> vtkspheremapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();

    vtkspheremapper->SetInputConnection(vtksphere->GetOutputPort());

    parent._multiblock->SetBlock(sphere.number(),
                                 vtkspheremapper->GetInput());
  }

};

struct SiconosVTKOutput::_DataSetUpdater : public SiconosVisitor
{
  SiconosVTKOutput& parent;

  _DataSetUpdater(SiconosVTKOutput& p) : parent(p) {};

  void visit(const SphereNEDS& sphere)
  {
    vtkSphereSource* vtksphere =
      static_cast<vtkSphereSource*>(parent._sources[sphere.number()]);

    vtksphere->SetCenter(sphere.q()->getValue(0),
                         sphere.q()->getValue(1),
                         sphere.q()->getValue(2));

    vtksphere->SetRadius(sphere.getRadius());
  }
};


SiconosVTKOutput::SiconosVTKOutput(SP::Model model, std::string filename) :
  SiconosOutput(model)
{
  _multiblock = vtkMultiBlockDataSet::New();
  _writer = vtkXMLMultiBlockDataWriter::New();
  _writer->SetFileName(filename.c_str());

  _DataSetMaker dataSetMaker(*this);

  SP::DynamicalSystemsGraph dsg = _model->
                                  nonSmoothDynamicalSystem()->topology()->dSG(0);

  DynamicalSystemsGraph::VIterator vi, viend;
  for (tie(vi, viend) = dsg->vertices();
       vi != viend; ++vi)
  {
    dsg->bundle(*vi)->accept(dataSetMaker);
  }

};

SiconosVTKOutput::~SiconosVTKOutput()
{
  if (_multiblock)
  {
    _multiblock->Delete();
  }

  if (_writer)
  {
    _writer->Delete();
  }
};


void SiconosVTKOutput::update()
{

  _DataSetUpdater dataSetUpdater = _DataSetUpdater(*this);

  SP::DynamicalSystemsGraph dsg = _model->
                                  nonSmoothDynamicalSystem()->topology()->dSG(0);

  // parallel IO -> use dfs or bfs algorithm
  DynamicalSystemsGraph::VIterator vi, viend;
  for (tie(vi, viend) = dsg->vertices();
       vi != viend; ++vi)
  {
    dsg->bundle(*vi)->accept(dataSetUpdater);
  }
}

void SiconosVTKOutput::write()
{
  _writer->Write();
}
