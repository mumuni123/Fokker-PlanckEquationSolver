#ifndef NONUNIFORM_RECONSTRUCTION_H
#define NONUNIFORM_RECONSTRUCTION_H

namespace NonuniformMuscl {

struct FaceStates {
    double left;
    double right;
};

// Monotonized-central MUSCL reconstruction of cell averages on physical,
// possibly nonuniform coordinates.  The four values are the two cells
// adjacent to the face and one physical neighbour on each side.
FaceStates reconstruct_face(double q_im1, double q_i, double q_ip1,
                            double q_ip2, double s_im1, double s_i,
                            double s_ip1, double s_ip2, double s_face);

// A face owns one Riemann state and therefore one conservative flux.
double upwind_state(const FaceStates& states, double speed);

// Symmetric face state used by the energy-compatible high-order u-force
// candidate.  It is deliberately separate from the donor/upwind state used
// by the low-order FCT safety transport.
double centered_state(const FaceStates& states);

// Converts a reconstructed cylindrical cell average into the conservative
// u_parallel-face coefficient C_u.  Species::f stores M=fbar*dx*du*A_perp,
// so C_u=M/du=fbar*dx*A_perp.
double upar_face_coefficient(double fbar, double dx, double transverse_area);

double upar_face_flux(double acceleration, double coefficient);

} // namespace NonuniformMuscl

#endif
