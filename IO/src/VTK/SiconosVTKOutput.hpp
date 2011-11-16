

#ifndef SICONOSVTKOUTPUT_HPP
#define SICONOSVTKOUTPUT_HPP

#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <SiconosKernel.hpp>
#include "SiconosOutput.hpp"

class SiconosVTKOutput : public SiconosOutput
{
protected:

  std::vector<vtkAlgorithm*> _sources;
  vtkMultiBlockDataSet* _multiblock;
  vtkXMLMultiBlockDataWriter* _writer;

  struct _DataSetMaker;
  struct _DataSetUpdater;

  friend class SiconosVTKOutput::_DataSetMaker;
  friend class SiconosVTKOutput::_DataSetUpdater;

public:

  SiconosVTKOutput() {};

  SiconosVTKOutput(SP::Model model, std::string filename);

  virtual ~SiconosVTKOutput();

  virtual void write();

  void update();

};



#endif
