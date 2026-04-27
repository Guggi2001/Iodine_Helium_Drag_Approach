# **How the potentials work**



Time t=0:

&#x20;   - Position: sitting inside the droplet at r ≈ 30 Å from center

&#x20;   - Partner: the other I atom in my I₂ is 2.666 Å away

&#x20;   - Forces on me:

&#x20;       (a) droplet\_potential(depth) → \~0 force (I'm deep inside)

&#x20;       (b) morse\_X(R) → strong binding force pulling me to partner



Laser fires! (dissociation):

&#x20;   - I now have kinetic energy from the photon

&#x20;   - My partner and I fly apart

&#x20;   - Forces on me:

&#x20;       (a) droplet\_potential → still \~0 (still inside)

&#x20;       (b) morse\_X(R) → residual repulsive force (Morse at R > R\_e)

&#x20;       (c) hard-sphere collisions with He atoms slow me down



At some later time, the atom approaches the droplet surface:

&#x20;   - Forces on me:

&#x20;       (a) droplet\_potential → strong force (I'm crossing the surface!)

&#x20;       (b) morse\_X(R) → now R is large, force is negligible

&#x20;       (c) hard-sphere collisions thin out



If the laser ionizes me:

&#x20;   - I'm now I⁺, partner is also I⁺

&#x20;   - Forces on me:

&#x20;       (a) droplet\_potential (different binding for ions)

&#x20;       (b) morse\_I2plus\_state\_select(R, my\_state) → strong!

&#x20;       (c) Coulomb repulsion

&#x20;       (d) He collisions (now with a different cross section)

