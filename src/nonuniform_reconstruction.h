#ifndef NONUNIFORM_RECONSTRUCTION_H
#define NONUNIFORM_RECONSTRUCTION_H

namespace NonuniformMuscl {

struct FaceStates {
    double left;
    double right;
};

struct FrozenCenteredLinearization {
    // Coefficients of the frozen centered state with respect to
    // (q_im1, q_i, q_ip1, q_ip2).  Branch 0 uses the zero-slope generalized
    // derivative, so this object remains finite at extrema and ties.
    double coefficient[4];
    int branch_signature;
};

// Exact linear stencil of one selected Riemann state after the MC/minmod
// choices have been frozen at the supplied state.  `use_left_state` must use
// the same velocity-sign selection as the production flux.  This is read-only
// audit support; it never changes the production reconstruction.
FrozenCenteredLinearization frozen_upwind_linearization(
    double q_im1, double q_i, double q_ip1, double q_ip2, double s_im1,
    double s_i, double s_ip1, double s_ip2, double s_face,
    bool use_left_state);

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

// Read-only audit identity for the two MC/minmod slopes used by
// reconstruct_face.  The low two bits encode the left-cell slope branch and
// the next two bits encode the right-cell slope branch: 0=zero/sign change,
// 1=centered candidate, 2=left candidate, 3=right candidate.  It is used by
// frozen-Jacobian audits only and does not alter reconstruction.
int centered_state_branch_signature(double q_im1, double q_i, double q_ip1,
                                    double q_ip2, double s_im1, double s_i,
                                    double s_ip1, double s_ip2);

// Freezes the MC/minmod choices selected at the supplied state and returns
// the exact local linear stencil for centered_state(reconstruct_face(...)).
// This is an audit operator; reconstruction itself is unchanged.
FrozenCenteredLinearization frozen_centered_linearization(
    double q_im1, double q_i, double q_ip1, double q_ip2, double s_im1,
    double s_i, double s_ip1, double s_ip2, double s_face);

// Converts a reconstructed cylindrical cell average into the conservative
// u_parallel-face coefficient C_u.  Species::f stores M=fbar*dx*du*A_perp,
// so C_u=M/du=fbar*dx*A_perp.
double upar_face_coefficient(double fbar, double dx, double transverse_area);

double upar_face_flux(double acceleration, double coefficient);

} // namespace NonuniformMuscl

#endif
