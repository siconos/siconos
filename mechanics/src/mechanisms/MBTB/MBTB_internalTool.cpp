#include "MBTB_DATA.hpp"
#include "MBTB_PYTHON_API.hpp"
#include "SiconosKernel.hpp"
#include "CADMBTB_API.hpp"
#include "ace.h"
#include "MBTB_TimeSteppingProj.hpp"
#include "MBTB_TimeSteppingCombinedProj.hpp"
#include "NewtonEulerJointR.hpp"

void _MBTB_updateContactFromDS()
{
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
#ifdef PRINT_FORCE_CONTACTS
    printf("....> contact force of %s :",sContacts[numC]->contactName());
    sContacts[numC]->relation()->contactForce()->display();
#endif
    int index1=sContacts[numC]->_indexBody1;
    int index2=sContacts[numC]->_indexBody2;

    CADMBTB_moveModelFromModel(sContacts[numC]->_indexCAD1,index1);
    if(index2!=-1)
      CADMBTB_moveModelFromModel(sContacts[numC]->_indexCAD2,index2);

#ifdef DRAW_CONTACT_FORCES

#endif
    CADMBTB_moveGraphicalModelFromModel(sContacts[numC]->_indexCAD1,index1);
    if(index2!=-1)
      CADMBTB_moveGraphicalModelFromModel(sContacts[numC]->_indexCAD2,index2);
  }
}




void _MBTB_updateContactFromDS(int numDS)
{
  for( int numC=0; numC < (int)sNbOfContacts; numC++)
  {
    if((int)(sContacts[numC]->_indexBody1) == numDS)
    {
      CADMBTB_moveModelFromModel(sContacts[numC]->_indexCAD1,numDS);
      CADMBTB_moveGraphicalModelFromModel(sContacts[numC]->_indexCAD1,numDS);
    }
    if(sContacts[numC]->_indexBody2 == numDS)
    {
      CADMBTB_moveModelFromModel(sContacts[numC]->_indexCAD2,numDS);
      CADMBTB_moveGraphicalModelFromModel(sContacts[numC]->_indexCAD2,numDS);
    }
  }
}

void _MBTB_DRAW_STEP()
{
  /*delete previous*/
  if(sDrawMode)
  {
    for(unsigned int nC=0; nC<3*sNbOfContacts; nC++)
    {
       CADMBTB_buildLineArtefactLine(nC,0,0,0,0,0,0);
    }

  }
  if(sDrawMode & MBTB_FACE_NORMAL1)
  {
    for(unsigned int nC=0; nC<sNbOfContacts; nC++)
    {
      double x1,x2,y1,y2,z1,z2;//,nx,ny,nz,MinDist;
      //int index1=sContacts[nC]->_indexBody1;
      //int index2=sContacts[nC]->_indexBody2;
      //int normalFromFace1;

      x2=sContacts[nC]->relation()->pc2()->getValue(0);
      y2=sContacts[nC]->relation()->pc2()->getValue(1);
      z2=sContacts[nC]->relation()->pc2()->getValue(2);
      x1=x2+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(0);
      y1=y2+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(1);
      z1=z2+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(2);
      /*  CADMBTB_getMinDistance(nC,index1,index2,
          x1,y1,z1,
          x2,y2,z2,
          nx,ny,nz,
          normalFromFace1,
          MinDist);
          printf ("second point of contact : x2=%lf ,y2=%lf,z2=%lf\n",x2,y2,z2);*/
      //  CADMBTB_buildOrientedLineArtefactLine(nC+sNbOfContacts,&x2,&y2,&z2,&x1,&y1,&z1);
    }
  }

  if(sDrawMode & MBTB_ARTEFACT_P1P2)
  {
    for(unsigned int nC=0; nC<sNbOfContacts; nC++)
    {
      double x1,x2,y1,y2,z1,z2;
      x1=sContacts[nC]->relation()->pc1()->getValue(0);
      y1=sContacts[nC]->relation()->pc1()->getValue(1);
      z1=sContacts[nC]->relation()->pc1()->getValue(2);
      x2=sContacts[nC]->relation()->pc2()->getValue(0);
      y2=sContacts[nC]->relation()->pc2()->getValue(1);
      z2=sContacts[nC]->relation()->pc2()->getValue(2);
      // printf ("second point of contact : x2=%lf ,y2=%lf,z2=%lf\n",x2,y2,z2);
      CADMBTB_buildLineArtefactLine(nC,&x1,&y1,&z1,&x2,&y2,&z2);

    }
  }


  if(sDrawMode & MBTB_ARTEFACT_NORMAL)
  {
    for(unsigned int nC=0; nC<sNbOfContacts; nC++)
    {
      double x1,x2,y1,y2,z1,z2;
      x1=sContacts[nC]->relation()->pc1()->getValue(0);
      y1=sContacts[nC]->relation()->pc1()->getValue(1);
      z1=sContacts[nC]->relation()->pc1()->getValue(2);
      x2=x1+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(0);
      y2=y1+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(1);
      z2=z1+1.1*sArtefactLength*sContacts[nC]->relation()->nc()->getValue(2);
      CADMBTB_buildOrientedLineArtefactLine(nC+sNbOfContacts,&x1,&y1,&z1,&x2,&y2,&z2);

    }
  }

  if(sDrawMode & MBTB_ARTEFACT_REACTION)
  {
    // printf("MBTB_DRAW_STEP REACTION\n");
    SP::Topology topo = sSimu->nonSmoothDynamicalSystem()->topology();
    double h = sSimu->timeStep();
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);


    double FMax=0;
    double aux,normF;
    int nbUR=0;
    int nbR=0;
    InteractionsGraph::VIterator ui1, ui1end;
    boost::tie(ui1,ui1end) = indexSet1->vertices();

    if(sNominalForce > 1e-12)
      FMax=sNominalForce;
    else
      for(;
          ui1 != ui1end; ++ui1)
      {
        nbUR++;
        SP::Interaction inter1 = indexSet1->bundle(*ui1);
        SP::Relation R=inter1->relation();
        for(unsigned int nC=0; nC<sNbOfContacts; nC++)
        {
          if(sContacts[nC]->relation() == R)
          {
            nbR++;
            SP::SiconosVector F=sContacts[nC]->relation()->contactForce();
            aux = sqrt(F->getValue(0)*F->getValue(0)+F->getValue(1)*F->getValue(1)+F->getValue(2)*F->getValue(2))/h;
            if(aux>FMax)
              FMax=aux;
          }
        }
      }
    //    printf("MBTB_DRAW_STEP REACTION nb Ur=%i, nb Rcontact=%i.\n",nbUR,nbR);

    // printf("_MBTB_DRAW_STEP FMAX=%e \n",FMax);

    boost::tie(ui1,ui1end) = indexSet1->vertices();

    for(;
        ui1 != ui1end; ++ui1)
    {
      SP::Interaction inter1 = indexSet1->bundle(*ui1);
      SP::Relation R=inter1->relation();
      for(unsigned int nC=0; nC<sNbOfContacts; nC++)
      {
        if(sContacts[nC]->relation() == R)
        {
          SP::SiconosVector F=sContacts[nC]->relation()->contactForce();
          normF=sqrt(F->getValue(0)*F->getValue(0)+F->getValue(1)*F->getValue(1)+F->getValue(2)*F->getValue(2));
          aux = normF/h;
          aux = aux/FMax;
          if(aux>1)
            aux=log(aux)+1.0;
          else
            aux=1.0/(-log(aux)+1.0);
          if(aux>5)
            aux=5;

          //printf("_MBTB_DRAW_STEP aux=%e \n",aux);
          if(aux>sArtefactThreshold)
          {
            SP::SiconosVector F=sContacts[nC]->relation()->contactForce();
            double x1,x2,y1,y2,z1,z2;
            x1=sContacts[nC]->relation()->pc1()->getValue(0);
            y1=sContacts[nC]->relation()->pc1()->getValue(1);
            z1=sContacts[nC]->relation()->pc1()->getValue(2);

            x2=x1+aux*(sArtefactLength/normF)*F->getValue(0);
            y2=y1+aux*(sArtefactLength/normF)*F->getValue(1);
            z2=z1+aux*(sArtefactLength/normF)*F->getValue(2);
            double radius =0.03*aux*sArtefactLength;
            CADMBTB_buildCylinderArtefactLine(nC+2*sNbOfContacts,&x1,&y1,&z1,&x2,&y2,&z2,&radius);
          }
        }
      }
    }
  }
  CADMBTB_updateGraphic();
}
void _MBTB_STEP()
{
  MBTB_updateDSFromSiconos();
  _MBTB_updateContactFromDS();
  ACE_times[ACE_TIMER_SICONOS].start();
  sSimu->setNewtonTolerance(sDParams[2]);
  sSimu->setNewtonMaxIteration(sDParams[3]);
  sSimu->advanceToEvent();

  double * dd =   sSimu->oneStepNSProblem(0)->numericsSolverOptions()->dparam;
  int * ii =   sSimu->oneStepNSProblem(0)->numericsSolverOptions()->iparam;

  std::cout<< "     OSNS reached accuracy ="<< dd[2] << " < " << dd [0] <<  std::endl;
  std::cout<< "     OSNS nb iterations ="<< ii[3] << " < " << ii [0] <<  std::endl;
  std::cout<< "     Number of Newton iterations = " << sSimu->getNewtonNbIterations() <<std::endl;

  Type::Siconos  simuType;
  simuType = Type::value(*sSimu);
  if (simuType == Type::TimeStepping)
  {
    
  }
  else if (simuType == Type::TimeSteppingDirectProjection)
  {
  }
  else if (simuType == Type::TimeSteppingCombinedProjection)
  {
    std::cout<< "     Number of projection iterations = " <<  (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->nbProjectionIteration() <<std::endl;
    std::cout<< "     Number of cumulated Newton iterations = " <<  (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->cumulatedNewtonNbIterations() <<std::endl;
    std::cout<< "     Number of set  iterations = " <<  (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->nbIndexSetsIteration() <<std::endl;
    std::cout<< "     Max violation unilateral = " <<  (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->maxViolationUnilateral()  <<std::endl;
    std::cout<< "     Max violation equality = " <<  (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->maxViolationEquality() <<std::endl;
  }
  
  //sSimu->oneStepNSProblem(0)->display();
  ACE_times[ACE_TIMER_SICONOS].stop();


  if(sTimerCmp%sFreqGraphic==0)
  {
    _MBTB_DRAW_STEP();
  }
  sSimu->nextStep();
}



void _MBTB_displayStep()
{
  //fprintf(fp,"_MBTB_printStep ");
  //printf("%d\t",sTimerCmp);

  // Bodies display output
  if (sDisplayStepBodies)
  {
    printf("STEP Number = %d\t",sTimerCmp);
    for(unsigned int numDS =0; numDS<sNbOfBodies; numDS++)
    {
      printf("Body number %i\n", numDS);
      printf("Position of body %i\n", numDS);
      for(unsigned int ii=0; ii<sDS[numDS]->q()->size(); ii++)
      {
        printf("%e",sDS[numDS]->q()->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Velocity of body %i\n", numDS);
      for(unsigned int ii=0; ii<sDS[numDS]->twist()->size(); ii++)
      {
        printf("%e",sDS[numDS]->twist()->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Kinetic Energy of body %i\n", numDS);
      /*Ec of the DS*/
      //   printf("MBTB Ec computattiom masse matrix:\n");
      //    (sDS[numDS]->M())->display();
      //    (sDS[numDS]->twist())->display();
      SiconosVector res(6);
      prod(*(sDS[numDS]->mass()),*(sDS[numDS]->twist()),res);
      double ec=0.0;
      for(int i=0; i<6; i++)
        ec+=res.getValue(i)*sDS[numDS]->twist()->getValue(i);
      printf("%e\t",ec*0.5);
      printf("\n");
    }
  }


  // Joints display output
  if (sDisplayStepJoints)
  {
    printf("STEP Number = %d\t",sTimerCmp);
    for(unsigned int numJ=0; numJ<sNbOfJoints; numJ++)
    {
      printf("Joint number %i\n", numJ);
      printf("interactionjointR->display  %i\n", numJ);
      sJointRelations[numJ]->_interaction->display();
      printf("\n");printf("Forces in Joint  %i\n", numJ);
      for(int ii=0; ii<3; ii++)
      {
        printf("%e",sJointRelations[numJ]->_jointR->contactForce()->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Moments in Joint  %i\n", numJ);

      SP::SiconosVector vaux(new SiconosVector(3));
      for(int ii=3; ii<6; ii++)
      {
        vaux->setValue(ii-3,sJointRelations[numJ]->_jointR->contactForce()->getValue(ii));
        printf("%e",vaux->getValue(ii-3));
        printf("\t");
      }
      /*convert momentum in abs frame*/
      printf("\n");

      printf("Moments in Joint %i in absolute frame \n", numJ);
      changeFrameBodyToAbs(sJointRelations[numJ]->_ds1->q(),vaux);
      for(int ii=0; ii<3; ii++)
      {
        printf("%e",vaux->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Equivalent Forces to moments in Joint %i in absolute frame \n", numJ);
      if(sJointRelations[numJ]->_G0C1)
      {
        sJointRelations[numJ]->computeEquivalentForces();
        for(int ii=0; ii <6; ii++)
        {
          printf("%e",sJointRelations[numJ]->_F->getValue(ii));
          printf("\t");
        }
      }
      else
      {
        printf("N/A\t");
      }
     printf("\n");

    }
  }
  // Contacts display output
  if (sDisplayStepContacts)
  {
    printf("STEP Number = %d\t",sTimerCmp);
    for(unsigned int numC=0; numC<sNbOfContacts; numC++)
    {
      printf("Contact number %i\n", numC);
      printf("Contact forces in contact  %i\n", numC);
      for(int ii=0; ii<3; ii++)
      {
        printf("%e",sContacts[numC]->relation()->contactForce()->getValue(ii));
        printf("\t");
      }


      SP::SiconosVector vaux(new SiconosVector(3));
      for(int ii=3; ii<6; ii++)
      {
        vaux->setValue(ii-3,sContacts[numC]->relation()->contactForce()->getValue(ii));
      }
      /*convert momentum in abs frame*/
      changeFrameBodyToAbs(sDS[sContacts[numC]->_indexBody1]->q(),vaux);
      printf("\n");
      printf("Moments of contact forces in contact  %i in absolute frame \n", numC);
      for(int ii=0; ii<3; ii++)
      {
        printf("%e",vaux->getValue(ii));
        printf("\t");
      }

      printf("\n");

      SP::InteractionsGraph indexSet1=myModel->nonSmoothDynamicalSystem()->topology()->indexSet(1);
      InteractionsGraph::VIterator ui1, ui1end, v1next;
      boost::tie(ui1,ui1end) = indexSet1->vertices();
      //Remove interactions from the indexSet1
      int find = 0;
      for(v1next = ui1 ;
          ui1 != ui1end; ui1 = v1next)
      {
        ++v1next;
        SP::Interaction urI = indexSet1->bundle(*ui1);
        if(&(*urI) == &(*(sInterContacts[numC])))
        {
          find =1;
        }
      }

      printf("Contact status %d\n",find);
      printf("Coordinates of first contact point for contact %i",numC);

      printf("%e\t%e\t%e\t",
             sContacts[numC]->relation()->pc1()->getValue(0),
             sContacts[numC]->relation()->pc1()->getValue(1),
             sContacts[numC]->relation()->pc1()->getValue(2));
      printf("\n");
      if (sContacts[numC]->relation()->pc2())
      {
        printf("Coordinates of second contact point for contact %i",numC);
        printf("%e\t%e\t%e\t",
               sContacts[numC]->relation()->pc2()->getValue(0),
               sContacts[numC]->relation()->pc2()->getValue(1),
               sContacts[numC]->relation()->pc2()->getValue(2));
        printf("\n");
      }
      printf("Gap  for contact %i\t =\t ",numC);
      printf("%e\n",
	     sContacts[numC]->interaction()->y(0)->getValue(0));
       printf("vitess  for contact %i\t =\t ",numC);
      printf("%e\n",
	     sContacts[numC]->interaction()->y(1)->getValue(0));

    }
  }
  // printf("\n");
}
FILE* _MBTB_open(std::string filename, std::string args)
{
  FILE* fp = fopen(filename.c_str(), args.c_str());
  return fp;
}

void _MBTB_close(FILE* fp)
{
  fclose(fp);
}


void _MBTB_printStep(FILE *fp)
{
  //fprintf(fp,"_MBTB_printStep ");
  fprintf(fp,"%d\t",sTimerCmp);
  for(unsigned int numDS =0; numDS<sNbOfBodies; numDS++)
  {
    for(unsigned int ii=0; ii<sDS[numDS]->q()->size(); ii++)
    {
      fprintf(fp,"%e",sDS[numDS]->q()->getValue(ii));
      fprintf(fp,"\t");
    }
    for(unsigned int ii=0; ii<sDS[numDS]->twist()->size(); ii++)
    {
      fprintf(fp,"%e",sDS[numDS]->twist()->getValue(ii));
      fprintf(fp,"\t");
    }
    /*Ec of the DS*/
    //   printf("MBTB Ec computattiom masse matrix:\n");
    //    (sDS[numDS]->M())->display();
    //    (sDS[numDS]->twist())->display();
    SiconosVector res(6);
    prod(*(sDS[numDS]->mass()),*(sDS[numDS]->twist()),res);
    double ec=0.0;
    for(int i=0; i<6; i++)
      ec+=res.getValue(i)*sDS[numDS]->twist()->getValue(i);
    fprintf(fp,"%e\t",ec*0.5);
    //    printf("ec=%e",ec*0.5);
  }
  for(unsigned int numJ=0; numJ<sNbOfJoints; numJ++)
  {
    for(unsigned int ii=0; ii<3; ii++)
    {
      fprintf(fp,"%e",sJointRelations[numJ]->_jointR->contactForce()->getValue(ii));
      fprintf(fp,"\t");
    }
    SP::SiconosVector vaux(new SiconosVector(3));
    for(unsigned int ii=3; ii<6; ii++)
    {
      vaux->setValue(ii-3,sJointRelations[numJ]->_jointR->contactForce()->getValue(ii));
    }
    /*convert momentum in abs frame*/
    changeFrameBodyToAbs(sJointRelations[numJ]->_ds1->q(),vaux);
    for(int ii=0; ii<3; ii++)
    {
      fprintf(fp,"%e",vaux->getValue(ii));
      fprintf(fp,"\t");
    }


    if(sJointRelations[numJ]->_G0C1)
    {
      sJointRelations[numJ]->computeEquivalentForces();
      for(int ii=0; ii <6; ii++)
      {
        fprintf(fp,"%e",sJointRelations[numJ]->_F->getValue(ii));
        fprintf(fp,"\t");
      }
    }

  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    for(unsigned int ii=0; ii<3; ii++)
    {
      fprintf(fp,"%e",sContacts[numC]->relation()->contactForce()->getValue(ii));
      fprintf(fp,"\t");
    }
    SP::SiconosVector vaux(new SiconosVector(3));
    for(unsigned int ii=3; ii<6; ii++)
    {
      vaux->setValue(ii-3,sContacts[numC]->relation()->contactForce()->getValue(ii));
    }
    /*convert momentum in abs frame*/
    changeFrameBodyToAbs(sDS[sContacts[numC]->_indexBody1]->q(),vaux);
    for(int ii=0; ii<3; ii++)
    {
      fprintf(fp,"%e",vaux->getValue(ii));
      fprintf(fp,"\t");
    }
  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {

    SP::InteractionsGraph indexSet1=myModel->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    InteractionsGraph::VIterator ui1, ui1end, v1next;
    boost::tie(ui1,ui1end) = indexSet1->vertices();
    //Remove interactions from the indexSet1
    int find = 0;
    for(v1next = ui1 ;
        ui1 != ui1end; ui1 = v1next)
    {
      ++v1next;
      SP::Interaction urI = indexSet1->bundle(*ui1);
      if(&(*urI) == &(*(sInterContacts[numC])))
      {
        find =1;
      }
    }
    fprintf(fp,"%d\t",find);

  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    fprintf(fp,"%e\t%e\t%e\t",
            sContacts[numC]->relation()->pc1()->getValue(0),
            sContacts[numC]->relation()->pc1()->getValue(1),
            sContacts[numC]->relation()->pc1()->getValue(2));

  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    unsigned int sizeY = sContacts[numC]->interaction()->y(0)->size();
    if (sizeY==1)
    {
      fprintf(fp,"%e\t0.\t0.\t",
              sContacts[numC]->interaction()->y(0)->getValue(0));
    }
    else
    {
    fprintf(fp,"%e\t%e\t%e\t",
            sContacts[numC]->interaction()->y(0)->getValue(0),
            sContacts[numC]->interaction()->y(0)->getValue(1),
            sContacts[numC]->interaction()->y(0)->getValue(2));
    }
  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    unsigned int sizeY = sContacts[numC]->interaction()->y(0)->size();
    if (sizeY==1)
    {
      fprintf(fp,"%e\t0.\t0.\t",
            sContacts[numC]->interaction()->y(1)->getValue(0));
    }
    else
    {
      fprintf(fp,"%e\t%e\t%e\t",
            sContacts[numC]->interaction()->y(1)->getValue(0),
            sContacts[numC]->interaction()->y(1)->getValue(1),
            sContacts[numC]->interaction()->y(1)->getValue(2));
    }
  }

 fprintf(fp,"\n");
}
void _MBTB_printHeader(FILE *fp)
{
  unsigned int cmp=1;
  //fprintf(fp,"_MBTB_printStep\t");
  fprintf(fp,"stepNum1\t");
  cmp++;
  for(unsigned int numDS =0; numDS<sNbOfBodies; numDS++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+7-1; ++icmp)
    {
      fprintf(fp,"position%d_ds_%d\t",icmp,numDS);
    }
    cmp+=7;
    for(unsigned int icmp=cmp; icmp<=cmp+6-1; ++icmp)
    {
      fprintf(fp,"velocity%d_ds_%d\t",icmp,numDS);
    }
    cmp+=6;
    fprintf(fp,"EC%d_ds_%d\t",cmp,numDS);
    cmp++;
  }
  for(unsigned int numJ=0; numJ<sNbOfJoints; numJ++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+6-1; ++icmp)
    {
      fprintf(fp,"jointF%d_%d\t",icmp,numJ);
    }
    cmp+=6;
    if(sJointRelations[numJ]->_G0C1)
    {
      for(unsigned int icmp=cmp; icmp<=cmp+6-1; ++icmp)
      {
        fprintf(fp,"jointEquiF%d_%d\t",icmp,numJ);
      }
      cmp+=6;
    }
  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+6-1; ++icmp)
    {
      fprintf(fp,"ContactForce_%s_%d,%d\t",sContacts[numC]->contactName(),icmp,numC);
    }
    cmp+=6;
  }
  /*  for(int numC=0; numC<sNbOfContacts; numC++)
  {
    fprintf(fp,"MomentsOfContactForce_%s_%d_%d\t",sContacts[numC]->contactName(),cmp,cmp+6-1,numC);
    cmp+=6;
    }*/
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    fprintf(fp,"ContactState_%s_%d\t",sContacts[numC]->contactName(),cmp);
    cmp++;
  }
  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+2; ++icmp)
    {
      fprintf(fp,"ContactPoint_%s_%d\t",sContacts[numC]->contactName(),icmp);
    }
    cmp+=3;
  }
 for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+2; ++icmp)
    {
      fprintf(fp,"ContactGap_%s_%d\t",sContacts[numC]->contactName(),icmp);
    }
    cmp+=3;
  }

 for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {
    for(unsigned int icmp=cmp; icmp<=cmp+2; ++icmp)
    {
      fprintf(fp,"ContactVitess_%s_%d\t",sContacts[numC]->contactName(),icmp);
    }
    cmp+=3;
  }
  fprintf(fp,"\n");
  _MBTB_printStep(fp);
}
