#include "CableCollisionManager.h"
#include "Cable2d3DR.hpp"
//#include "NonSmoothLaw.h"

CableCollisionManager::CableCollisionManager(const std::shared_ptr<CableDS> a_model,
                        const vector<std::shared_ptr<Support>> &a_supports)
{
	m_model = a_model;
	m_supports = a_supports;
	// boucle sur les supports (pulleys, piles)
	// récupère les positions de contact (fonction du maillage)
	// ??

    // ==> à mettre dans les suppports dépend si piles ou pulley
	unsigned int dim = 2; //( friction 2D)
    //shared_ptr<NonSmoothLaw> nslaw1 = make_shared<NonSmoothLaw>(e, mu_S, eT, dim);
    //shared_ptr<NonSmoothLaw> nslaw2 = make_shared<NonSmoothLaw>(e, mu_p, eT, dim);

    
}

CableCollisionManager::~CableCollisionManager() {}


/** Called by Simulation after updating positions prior to starting
 * the Newton loop. */
void CableCollisionManager::updateInteractions(std::shared_ptr < Simulation> simulation) 
{
  SiconosVector &q = *(m_model->q());
  size_t nb = q.size();

  unsigned int node_idx=0;
  for (size_t i=0; i<nb; i+=3) { // boucle par 3 pour récupérer x,y,z
    ContactPoints *contact = nullptr;// findContactPoint(node_idx);
    for (auto &s : m_supports) {
      if (s->isContact(Point(q.getValue(i), q.getValue(i+1), q.getValue(i+2)), m_tolContact)) { // test si le point x,y,z est en contact avec lel support
        // récupérer pc2 (projection du point sur l'obstacle), normal, tangent
        //s.getContact(x, y, z)
        
        if (contact) {
            // r = contact->getRelation
            // r->updateContactPoint(pc1, pc2, normal, tangent);
        }       
        else {
          // create relation
          std::shared_ptr<Cable2d3DR> r = make_shared<Cable2d3DR>(node_idx);
          //r->updateContactPoint(pc1, pc2, normal, tangent);

          // create interaction
          // std::shared_ptr<Interaction> inter = std::make_shared<Interaction>(_nslaw, r);
          // link the interaction and the dynamical system
          // simulation->link(m_model, inter);    
        }        
      }
      else {
        // pas en contact
        // si existe un contact -> remove        
      }
    }
  }
}
