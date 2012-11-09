

#ifndef SICONOSVTKOUTPUT_HPP
#define SICONOSVTKOUTPUT_HPP

#include <vtkMultiBlockDataSet.h>
#include <vtkTemporalDataSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLMultiBlockDataWriter.h>
#include <SiconosKernel.hpp>
#include "SiconosOutput.hpp"

class SiconosVTKOutput : public SiconosOutput
{
protected:

  vtkMultiBlockDataSet* _multiblock;
  vtkTemporalDataSet* _temporal_data;
  vtkXMLMultiBlockDataWriter* _writer;

  struct _DataSetMaker;

  friend class SiconosVTKOutput::_DataSetMaker;

public:

  SiconosVTKOutput() {};

  SiconosVTKOutput(SP::Model model, std::string filename);

  virtual ~SiconosVTKOutput();

  void write();

};



#endif
