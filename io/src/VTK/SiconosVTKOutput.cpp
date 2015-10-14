#include "SiconosVTKOutput.hpp"
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>

#include <SphereNEDS.hpp>


struct SiconosVTKOutput::_DataSetMaker : public SiconosVisitor
{
  SiconosVTKOutput& parent;
  vtkSmartPointer<vtkStringArray> shape_info;
  vtkSmartPointer<vtkDoubleArray> attributes;
  vtkSmartPointer<vtkDoubleArray> translation;
  vtkSmartPointer<vtkDoubleArray> orientation;

  _DataSetMaker(SiconosVTKOutput& p) : parent(p)
  {
    shape = vtkSmartPointer<vtkStringArray>::New();
    shape->SetName('shape');
    attributes = vtkSmartPointer<vtkDoubleArray>::New();
    attributes->SetName('attributes');
    attributes->SetNumberOfComponents(3);
    translation = vtkSmartPointer<vtkDoubleArray>::New();
    translation->SetName('translation');
    translation->SetNumberOfComponents(3);
    orientation = vtkSmartPointer<vtkDoubleArray>::New();
    orientation->SetName('orientation');
    orientation->SetNumberOfComponents(4);
  };

  void visit(const SphereNEDS& sphere)
  {

    shape_info->InsertNextValue('S');
    attributes->InsertNextTuple3(sphere->radius(), 0., 0.);

    translation->InsertNextTuple3(sphere.q()->getValue(0),
                                  sphere.q()->getValue(1),
                                  sphere.q()->getValue(2));

    orientation->InsertNextTuple4(sphere.q()->getValue(0),
                                  sphere.q()->getValue(1),
                                  sphere.q()->getValue(2),
                                  sphere.q()->getValue(3),
                                  sphere.q()->getValue(4));


  }

  void setTime(double time)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

    grid->GetPointData()->AddArray(shape);
    grid->GetPointData()->AddArray(attributes);
    grid->GetPointData()->AddArray(translation);
    grid->GetPointData()->AddArray(rotation);

    parent._temporal_data.SetTimeStep(time, grid);
  }

};


SiconosVTKOutput::SiconosVTKOutput(SP::Model model, std::string filename) :
  SiconosOutput(model)
{
  _multiblock = vtkMultiBlockDataSet::New();
  _writer = vtkXMLMultiBlockDataWriter::New();
  _writer->SetFileName(filename.c_str());

  _DataSetMaker dataSetMaker(*this);
  _DataSetUpdater dataSetUpdater = _DataSetUpdater(*this);


  SP::DynamicalSystemsGraph dsg = _model->
                                  nonSmoothDynamicalSystem()->topology()->dSG(0);

  DynamicalSystemsGraph::VIterator vi, viend;
  for (tie(vi, viend) = dsg->vertices();
       vi != viend; ++vi)
  {
    dsg->bundle(*vi)->accept(dataSetMaker);
  }

  dataSetMaker.setTime(_model->currentTime());

};

SiconosVTKOutput::~SiconosVTKOutput()
{
  if (_multiblock)
  {
    _multiblock->Delete();
  }

  if (_temporal_data)
  {
    _temporal_data->Delete();
  }


  if (_writer)
  {
    _writer->Delete();
  }
};

void SiconosVTKOutput::write()
{
  _writer->Write();
}
