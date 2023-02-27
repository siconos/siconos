#include "CableCollisionManager.h"

#include "Cable2d3DR.hpp"
#include "Simulation.hpp"
// #include "NonSmoothLaw.h"

CableCollisionManager::CableCollisionManager(
    const std::shared_ptr<CableDS> a_model, const vector<std::shared_ptr<Support>> &a_supports,
    double a_tolContact)
{
  m_model = a_model;
  m_supports = a_supports;
  m_tolContact = a_tolContact;

  // ==> à mettre dans les suppports dépend si piles ou pulley
  unsigned int dim = 2; //( friction 2D)
  // shared_ptr<NonSmoothLaw> nslaw1 = make_shared<NonSmoothLaw>(e, mu_S, eT, dim);
  // shared_ptr<NonSmoothLaw> nslaw2 = make_shared<NonSmoothLaw>(e, mu_p, eT, dim);
}

CableCollisionManager::~CableCollisionManager() {}

/** Called by Simulation after updating positions prior to starting
 * the Newton loop. */
void CableCollisionManager::updateInteractions(std::shared_ptr<Simulation> simulation)
{
  SiconosVector &q = *(m_model->q());
  size_t nb = q.size();
  
  unsigned int node_idx = 0;
  for (size_t i = 0; i < nb; i += 3, node_idx++) { // boucle par 3 pour récupérer x,y,z

    t_contacts::iterator contactItr = m_contacts.find(node_idx);
    std::shared_ptr<Interaction> contact = nullptr;
    if (contactItr != m_contacts.end()) {
      contact = contactItr->second;
    }
    auto pc1 = std::make_shared<SiconosVector>(3);
    pc1->setValue(0, q.getValue(i));
    pc1->setValue(1, q.getValue(i + 1));
    pc1->setValue(2, q.getValue(i + 2));
    
    for (auto &s : m_supports) {

      if (s->isContact(pc1, m_tolContact)) {
	// test si le point x,y,z est en contact avec le support
        // récupérer pc2 (projection du point sur l'obstacle), normal, tangent
        std::shared_ptr<SiconosVector> pc2 = s->pc2();
        std::shared_ptr<SiconosVector> normal = s->normal();
        std::shared_ptr<SiconosVector> tangent = s->tangent();

        if (contact) {
          std::shared_ptr<Cable2d3DR> relation =
              std::static_pointer_cast<Cable2d3DR>(contact->relation());
          relation->updateContactPoint(pc1, pc2, normal, tangent);
        }
        else {
          // create relation
          auto relation = std::make_shared<Cable2d3DR>(node_idx, pc2, normal, tangent);

          // create interaction
          auto inter = std::make_shared<Interaction>(s->nslaw(), relation);
          m_contacts[node_idx] = inter;

          // link the interaction and the dynamical system
          simulation->link(inter, m_model);
        }
      }
      else {
        // pas en contact
        if (contact) {
          // si existe un contact -> remove
          m_contacts.erase(node_idx);
        }
      }
    }
  }
  std::cout << "Fin UpInterac \n"; 
}
