This repository contains minor modifications to Bullet functions, which we implement by deriving new classes.

- Contact detection between spheres is not persistent for contact points. Whilst the contact manifold is preserved by Bullet, the `clearManifold();` function is called each time a collision detection call is made. This results in the creation of new contact points and their associated interactions.
