#include "tail_subcell_quadrature.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

void gauss_legendre_rule(int points, std::vector<double>& x,
                         std::vector<double>& w)
{
    // Fixed positive rules avoid a numerical root finder in this diagnostic
    // path.  Four points are the production audit default.
    if (points <= 2) {
        x.assign(2, 0.0); w.assign(2, 1.0);
        x[0] = -0.57735026918962576451; x[1] = -x[0];
        return;
    }
    x.assign(4, 0.0); w.assign(4, 0.0);
    x[0] = -0.86113631159405257522; x[1] = -0.33998104358485626480;
    x[2] = -x[1]; x[3] = -x[0];
    w[0] = 0.34785484513745385737; w[1] = 0.65214515486254614263;
    w[2] = w[1]; w[3] = w[0];
}

double kinetic_energy(double upar, double uperp)
{
    const double gamma = std::sqrt(1.0 + upar * upar + uperp * uperp);
    return Const::me * Const::c * Const::c * (gamma - 1.0);
}

double momentum_radius_squared(double energy)
{
    const double normalized = energy / (Const::me * Const::c * Const::c);
    return std::max(0.0, normalized * (normalized + 2.0));
}

void append_symmetric_root(std::vector<double>& points, double squared,
                           double lo, double hi)
{
    if (!(squared > 0.0)) return;
    const double root = std::sqrt(squared);
    if (-root > lo && -root < hi) points.push_back(-root);
    if (root > lo && root < hi) points.push_back(root);
}

double shell_section_integral(double upar_lo, double upar_hi,
                              double uperp_lo_sq, double uperp_hi_sq,
                              double radius_lo_sq, double radius_hi_sq)
{
    std::vector<double> points;
    points.push_back(upar_lo);
    points.push_back(upar_hi);
    if (upar_lo < 0.0 && upar_hi > 0.0) points.push_back(0.0);
    append_symmetric_root(points, radius_lo_sq - uperp_lo_sq,
                          upar_lo, upar_hi);
    append_symmetric_root(points, radius_lo_sq - uperp_hi_sq,
                          upar_lo, upar_hi);
    append_symmetric_root(points, radius_hi_sq - uperp_lo_sq,
                          upar_lo, upar_hi);
    append_symmetric_root(points, radius_hi_sq - uperp_hi_sq,
                          upar_lo, upar_hi);
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end(),
        [](double a, double b) {
            return std::fabs(a - b) <=
                32.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
        }), points.end());

    double integral = 0.0;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        const double a = points[i];
        const double b = points[i + 1];
        if (!(b > a)) continue;
        const double mid = 0.5 * (a + b);
        const double mid_sq = mid * mid;
        const bool upper_is_cell =
            uperp_hi_sq <= radius_hi_sq - mid_sq;
        const bool lower_is_cell =
            uperp_lo_sq >= radius_lo_sq - mid_sq;
        const double upper_a = upper_is_cell ? uperp_hi_sq : radius_hi_sq;
        const double upper_b = upper_is_cell ? 0.0 : -1.0;
        const double lower_a = lower_is_cell ? uperp_lo_sq : radius_lo_sq;
        const double lower_b = lower_is_cell ? 0.0 : -1.0;
        const double value_mid =
            (upper_a - lower_a) + (upper_b - lower_b) * mid_sq;
        if (!(value_mid > 0.0)) continue;
        const double dx = b - a;
        const double cubic = (b * b * b - a * a * a) / 3.0;
        integral += (upper_a - lower_a) * dx +
                    (upper_b - lower_b) * cubic;
    }
    return std::max(0.0, integral);
}

} // namespace

std::vector<TailSubcellNode> TailSubcellQuadrature::nodes(
    const CylindricalVelocityGrid& grid, int iv, int imu,
    int points_per_dimension)
{
    std::vector<TailSubcellNode> result;
    if (iv < 0 || imu < 0 ||
        iv + 1 >= static_cast<int>(grid.upar_faces.size()) ||
        imu + 1 >= static_cast<int>(grid.uperp_faces.size())) return result;
    std::vector<double> x, w;
    gauss_legendre_rule(points_per_dimension, x, w);
    const double up_lo = grid.upar_faces[static_cast<size_t>(iv)];
    const double up_hi = grid.upar_faces[static_cast<size_t>(iv) + 1];
    const double ut_lo = grid.uperp_faces[static_cast<size_t>(imu)];
    const double ut_hi = grid.uperp_faces[static_cast<size_t>(imu) + 1];
    const double up_mid = 0.5 * (up_lo + up_hi);
    const double up_half = 0.5 * (up_hi - up_lo);
    const double ut_mid = 0.5 * (ut_lo + ut_hi);
    const double ut_half = 0.5 * (ut_hi - ut_lo);
    double normalization = 0.0;
    for (size_t a = 0; a < x.size(); ++a) {
        for (size_t b = 0; b < x.size(); ++b) {
            const double uperp = ut_mid + ut_half * x[b];
            normalization += w[a] * w[b] * std::max(0.0, uperp);
        }
    }
    if (!(normalization > 0.0)) return result;
    result.reserve(x.size() * x.size());
    double fraction_sum = 0.0;
    for (size_t a = 0; a < x.size(); ++a) {
        for (size_t b = 0; b < x.size(); ++b) {
            TailSubcellNode node;
            node.upar = up_mid + up_half * x[a];
            node.uperp = ut_mid + ut_half * x[b];
            node.mass_fraction = w[a] * w[b] * std::max(0.0, node.uperp) /
                                 normalization;
            node.kinetic_energy = kinetic_energy(node.upar, node.uperp);
            fraction_sum += node.mass_fraction;
            result.push_back(node);
        }
    }
    // Make the partition of a cell mass exact in floating point while
    // preserving non-negative weights.
    if (!result.empty()) result.back().mass_fraction += 1.0 - fraction_sum;
    return result;
}

int TailSubcellQuadrature::energy_bin(const std::vector<double>& edges,
                                      double energy)
{
    if (edges.size() < 2 || energy < edges.front() ||
        energy > edges.back()) return -1;
    if (energy == edges.back()) return static_cast<int>(edges.size()) - 2;
    return static_cast<int>(std::upper_bound(edges.begin(), edges.end(), energy) -
                            edges.begin()) - 1;
}

std::vector<TailEnergyBinFraction>
TailSubcellQuadrature::energy_bin_fractions(
    const CylindricalVelocityGrid& grid, int iv, int imu,
    const std::vector<double>& energy_edges)
{
    std::vector<TailEnergyBinFraction> result;
    if (iv < 0 || imu < 0 || energy_edges.size() < 2 ||
        iv + 1 >= static_cast<int>(grid.upar_faces.size()) ||
        imu + 1 >= static_cast<int>(grid.uperp_faces.size())) return result;
    for (size_t i = 1; i < energy_edges.size(); ++i) {
        if (!(energy_edges[i] > energy_edges[i - 1]) ||
            !std::isfinite(energy_edges[i])) return result;
    }

    const double up_lo = grid.upar_faces[static_cast<size_t>(iv)];
    const double up_hi = grid.upar_faces[static_cast<size_t>(iv) + 1];
    const double ut_lo = grid.uperp_faces[static_cast<size_t>(imu)];
    const double ut_hi = grid.uperp_faces[static_cast<size_t>(imu) + 1];
    const double ut_lo_sq = ut_lo * ut_lo;
    const double ut_hi_sq = ut_hi * ut_hi;
    const double denominator = (up_hi - up_lo) * (ut_hi_sq - ut_lo_sq);
    if (!(denominator > 0.0)) return result;

    const double min_abs_up = up_lo <= 0.0 && up_hi >= 0.0 ? 0.0 :
        std::min(std::fabs(up_lo), std::fabs(up_hi));
    const double max_abs_up =
        std::max(std::fabs(up_lo), std::fabs(up_hi));
    const double min_energy = kinetic_energy(min_abs_up, ut_lo);
    const double max_energy = kinetic_energy(max_abs_up, ut_hi);
    std::vector<double>::const_iterator first_it =
        std::upper_bound(energy_edges.begin(), energy_edges.end(), min_energy);
    size_t first_bin = first_it == energy_edges.begin() ? 0 :
        static_cast<size_t>(first_it - energy_edges.begin() - 1);
    const size_t bin_count = energy_edges.size() - 1;
    if (first_bin >= bin_count) return result;

    double fraction_sum = 0.0;
    for (size_t bin = first_bin; bin < bin_count; ++bin) {
        if (energy_edges[bin] > max_energy) break;
        const double radius_lo_sq =
            momentum_radius_squared(energy_edges[bin]);
        const double radius_hi_sq =
            momentum_radius_squared(energy_edges[bin + 1]);
        const double overlap = shell_section_integral(
            up_lo, up_hi, ut_lo_sq, ut_hi_sq,
            radius_lo_sq, radius_hi_sq);
        if (!(overlap > 0.0)) continue;
        TailEnergyBinFraction entry;
        entry.bin = static_cast<int>(bin);
        entry.mass_fraction = overlap / denominator;
        fraction_sum += entry.mass_fraction;
        result.push_back(entry);
    }

    const bool fully_covered = energy_edges.front() <= min_energy &&
                               energy_edges.back() >= max_energy;
    if (fully_covered && !result.empty()) {
        result.back().mass_fraction += 1.0 - fraction_sum;
        if (result.back().mass_fraction < 0.0 &&
            result.back().mass_fraction > -1.0e-13)
            result.back().mass_fraction = 0.0;
    }
    return result;
}

void TailSubcellQuadrature::accumulate_cell(
    const CylindricalVelocityGrid& grid, int iv, int imu, double cell_mass,
    const std::vector<double>& energy_edges, double threshold_energy,
    TailSubcellSpectrum& spectrum, int points_per_dimension)
{
    if (!(cell_mass >= 0.0) || !std::isfinite(cell_mass)) return;
    const std::vector<TailSubcellNode> q =
        nodes(grid, iv, imu, points_per_dimension);
    if (q.empty()) return;
    spectrum.input_center_energy +=
        cell_mass * grid.kinetic_energy[idx2(iv, imu)];
    bool below = false;
    bool above = false;
    for (size_t n = 0; n < q.size(); ++n) {
        const double mass = cell_mass * q[n].mass_fraction;
        const int bin = energy_bin(energy_edges, q[n].kinetic_energy);
        spectrum.input_number += mass;
        spectrum.represented_number += mass;
        spectrum.represented_energy += mass * q[n].kinetic_energy;
        below = below || q[n].kinetic_energy < threshold_energy;
        above = above || q[n].kinetic_energy >= threshold_energy;
        if (bin >= 0 && bin < static_cast<int>(spectrum.number.size())) {
            spectrum.number[static_cast<size_t>(bin)] += mass;
            spectrum.energy[static_cast<size_t>(bin)] +=
                mass * q[n].kinetic_energy;
        }
    }
    if (below && above) {
        ++spectrum.straddling_cell_count;
        spectrum.straddling_cell_mass += cell_mass;
    }
    spectrum.number_residual = spectrum.represented_number - spectrum.input_number;
    // This is a representation delta, not a conservation error: it measures
    // the kinetic-energy change between point-mass-at-centre interpretation
    // and the cell-constant finite-volume interpretation.
    spectrum.energy_residual =
        spectrum.represented_energy - spectrum.input_center_energy;
}
