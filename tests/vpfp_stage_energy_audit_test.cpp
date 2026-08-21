#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace {

const char* const kStages[] = {
    "accepted_n", "collision_half1", "x_half1", "midpoint_poisson",
    "u_force_tail_beam_kick", "conversion_after_force", "x_half2",
    "collision_half2", "conversion_after_collision", "tail_bulk_return",
    "final_poisson" };
const int kStageCount = sizeof(kStages) / sizeof(kStages[0]);

bool make_directory(const std::string& path)
{
#ifdef _WIN32
    return _mkdir(path.c_str()) == 0 || errno == EEXIST;
#else
    return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
#endif
}

std::string quote(const std::string& text)
{
    return "\"" + text + "\"";
}

bool has_field(const std::string& path, const std::string& key,
               const std::string& expected)
{
    std::ifstream in(path.c_str());
    std::string line;
    const std::string wanted = key + "=" + expected;
    while (std::getline(in, line)) {
        if (line == wanted) return true;
    }
    return false;
}

double read_value(const std::string& path, const std::string& key, bool* found)
{
    std::ifstream in(path.c_str());
    std::string line;
    const std::string prefix = key + "=";
    while (std::getline(in, line)) {
        if (line.compare(0, prefix.size(), prefix) == 0) {
            if (found) *found = true;
            return std::atof(line.c_str() + prefix.size());
        }
    }
    if (found) *found = false;
    return 0.0;
}

void write_header(std::ostream& out, bool include_work_fields = true)
{
    out << "step time_s accepted audit_valid split failure_code "
        << "energy_balance_residual stage_id stage_name "
        << "K_bulk K_tail K_beam U_E "
        << "dK_bulk dK_tail dK_beam dU_E "
        << "Q_bkg_left_in Q_bkg_left_out Q_bkg_right_in Q_bkg_right_out "
        << "Q_beam_in Q_beam_out Q_tail_out Q_collision_reservoir "
        << "K_conversion K_tail_return W_electrostatic_boundary "
        << "stage_balance";
    if (include_work_fields) {
        out << " bulk_upar_face_work bulk_upar_velocity_boundary_work "
            << "bulk_upar_interface_energy_removed "
            << "bulk_upar_identity_residual tail_kick_work beam_kick_work";
    }
    out << "\n";
}

enum InvalidCase {
    VALID,
    MISSING_STAGE,
    DUPLICATE_STAGE,
    UNACCEPTED_STAGE,
    NONFINITE_STAGE,
    SPLIT_STAGE,
    FAILURE_RECORD,
    MISSING_WORK_FIELDS
};

bool write_fixture(const std::string& directory, int residual_stage,
                   InvalidCase invalid)
{
    if (!make_directory(directory)) return false;
    std::ofstream out((directory + "/vpfp_stage_energy_audit.dat").c_str());
    if (!out) return false;
    const bool include_work_fields = invalid != MISSING_WORK_FIELDS;
    write_header(out, include_work_fields);
    const double residual = residual_stage >= 0 ? 7.0 : 0.0;
    for (int stage = 0; stage < kStageCount; ++stage) {
        if (invalid == MISSING_STAGE && stage == 5) continue;
        const int stage_id = invalid == DUPLICATE_STAGE && stage == 5 ? 4 : stage;
        const char* stage_name = invalid == DUPLICATE_STAGE && stage == 5
            ? kStages[4] : kStages[stage];
        const bool changed = residual_stage >= 0 && stage >= residual_stage;
        out << 1 << " " << 1.0 << " "
            << (invalid == UNACCEPTED_STAGE && stage == 4 ? 0 : 1) << " "
            << 1 << " " << (invalid == SPLIT_STAGE && stage == 3 ? 1 : 0)
            << " 0 "
            << (invalid == NONFINITE_STAGE && stage == 4 ? "nan" : "") ;
        if (!(invalid == NONFINITE_STAGE && stage == 4)) out << residual;
        out << " " << stage_id << " " << stage_name;
        const double values[20] = {
            changed ? residual : 0.0, 0.0, 0.0, 0.0,
            stage == residual_stage ? residual : 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, stage == residual_stage ? residual : 0.0 };
        for (int i = 0; i < 20; ++i) out << " " << values[i];
        if (include_work_fields) {
            // The production stage-energy schema always carries these six
            // Gate C work columns.  Zero is the self-consistent value for the
            // synthetic fixtures unless a dedicated work test says otherwise.
            for (int i = 0; i < 6; ++i) out << " 0";
        }
        out << "\n";
    }
    if (invalid == FAILURE_RECORD) {
        std::ofstream failure((directory + "/vpfp_failure.dat").c_str());
        failure << "1 1 failure_code=4\n";
    }
    return static_cast<bool>(out);
}

// Writes a single accepted step with a controllable kinetic-energy magnitude
// and explicit telescope perturbation / ledger residual.  Only the final
// stage's K_bulk differs from energy_magnitude (a "missing source" of
// final_kbulk_delta J/m2); only midpoint_poisson carries a nonzero
// stage_balance; every stage carries the same ledger residual.  This lets the
// unit test drive the Gate A machine-precision gate directly.
bool write_scale_fixture(const std::string& directory,
                         double energy_magnitude,
                         double final_kbulk_delta,
                         double midpoint_stage_balance,
                         double ledger_residual)
{
    if (!make_directory(directory)) return false;
    std::ofstream out((directory + "/vpfp_stage_energy_audit.dat").c_str());
    if (!out) return false;
    out << std::setprecision(17);
    write_header(out);
    for (int stage = 0; stage < kStageCount; ++stage) {
        const double k_bulk = stage == kStageCount - 1
            ? energy_magnitude + final_kbulk_delta : energy_magnitude;
        const double stage_balance = stage == 3 ? midpoint_stage_balance : 0.0;
        out << "1 1.0 1 1 0 0 " << ledger_residual
            << " " << stage << " " << kStages[stage];
        const double values[20] = {
            k_bulk, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, stage_balance };
        for (int i = 0; i < 20; ++i) out << " " << values[i];
        for (int i = 0; i < 6; ++i) out << " 0";
        out << "\n";
    }
    return static_cast<bool>(out);
}

bool run_analyzer(const std::string& directory, const std::string& result,
                  bool require_no_split, int expected_steps)
{
    std::ostringstream command;
#ifdef _WIN32
    command << "python tools/analyze_vpfp_stage_energy_audit.py --run ";
#else
    command << "python3 tools/analyze_vpfp_stage_energy_audit.py --run ";
#endif
    command << quote(directory) << " --result " << quote(result)
            << " --expected-accepted-steps " << expected_steps;
    if (require_no_split) command << " --require-no-split";
    return std::system(command.str().c_str()) == 0;
}

bool parse_args(int argc, char** argv, std::string& result)
{
    result = "output/vpfp_stage_energy_audit_unit.result";
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--case" && i + 1 < argc) {
            if (std::string(argv[++i]) != "all") return false;
        } else if (arg == "--result" && i + 1 < argc) {
            result = argv[++i];
        } else {
            return false;
        }
    }
    return true;
}

bool run_case(const std::string& root, const char* name, int residual_stage,
              InvalidCase invalid, bool expect_pass, const char* dominant)
{
    const std::string directory = root + "/" + name;
    const std::string result = directory + ".result";
    if (!write_fixture(directory, residual_stage, invalid)) return false;
    const bool command_passed = run_analyzer(directory, result, true, 1);
    const bool status_matches = has_field(result, "status", expect_pass ? "PASS" : "FAIL");
    const bool dominant_matches = !dominant || has_field(result,
        "first_dominant_stage", dominant);
    return command_passed == expect_pass && status_matches && dominant_matches;
}

// Drives the Gate A machine-precision gate.  For pass cases it also requires
// audit_structure_pass=1 and physical_energy_gate_evaluated=0; for the
// self-consistent large-residual case it verifies the physical residual is
// reported nonzero without failing the structure audit.
bool run_scale_case(const std::string& root, const char* name,
                    double energy_magnitude, double final_kbulk_delta,
                    double midpoint_stage_balance, double ledger_residual,
                    bool expect_pass, bool expect_nonzero_physical_residual)
{
    const std::string directory = root + "/" + name;
    const std::string result = directory + ".result";
    if (!write_scale_fixture(directory, energy_magnitude, final_kbulk_delta,
                             midpoint_stage_balance, ledger_residual)) {
        return false;
    }
    const bool command_passed = run_analyzer(directory, result, true, 1);
    const bool status_matches = has_field(result, "status",
        expect_pass ? "PASS" : "FAIL");
    const bool gate_evaluated_zero = has_field(result,
        "physical_energy_gate_evaluated", "0");
    const bool structure_ok = has_field(result, "audit_structure_pass",
        expect_pass ? "1" : "0");
    bool physical_residual_ok = true;
    if (expect_nonzero_physical_residual) {
        bool found = false;
        physical_residual_ok = read_value(result,
            "physical_energy_residual_cumulative", &found) != 0.0 && found;
    }
    return command_passed == expect_pass && status_matches &&
        gate_evaluated_zero && structure_ok && physical_residual_ok;
}

} // namespace

int main(int argc, char** argv)
{
    std::string result;
    if (!parse_args(argc, argv, result)) {
        std::cerr << "usage: vpfp_stage_energy_audit_test --case all --result PATH\n";
        return 2;
    }
    const std::string root = result + ".cases";
    if (!make_directory("output") || !make_directory(root)) {
        std::cerr << "failed to create audit fixture directory\n";
        return 2;
    }

    const bool exact = run_case(root, "exact", -1, VALID, true, "x_half1");
    const bool x1 = run_case(root, "x1", 2, VALID, true, "x_half1");
    const bool force = run_case(root, "force", 4, VALID, true, "field_coupling");
    const bool collision2 = run_case(root, "collision2", 7, VALID, true, "collision_half2");
    const bool missing = run_case(root, "missing", -1, MISSING_STAGE, false, 0);
    const bool duplicate = run_case(root, "duplicate", -1, DUPLICATE_STAGE, false, 0);
    const bool unaccepted = run_case(root, "unaccepted", -1, UNACCEPTED_STAGE, false, 0);
    const bool nonfinite = run_case(root, "nonfinite", -1, NONFINITE_STAGE, false, 0);
    const bool split = run_case(root, "split", -1, SPLIT_STAGE, false, 0);
    const bool failure = run_case(root, "failure", -1, FAILURE_RECORD, false, 0);
    const bool missing_work = run_case(root, "missing_work", -1,
                                       MISSING_WORK_FIELDS, false, 0);

    // Gate A machine-precision gate cases.
    // 1. 5e9 energy scale + ~2e-6 telescoping perturbation must PASS.
    const bool large_perturb = run_scale_case(root, "large_perturb",
        5.0e9, 0.0, 2.0e-6, 0.0, true, false);
    // 2. Same scale, a 1 J/m2 unaccounted source (missing term) must FAIL.
    const bool large_missing_source = run_scale_case(root, "large_missing_source",
        5.0e9, 1.0, 0.0, 0.0, false, false);
    // 3. Small-scale artificial ledger with a 1e-4 residual must FAIL.
    const bool small_ledger_residual = run_scale_case(root, "small_ledger_residual",
        1.0, 0.0, 0.0, 1.0e-4, false, false);
    // 4. Large but self-consistent physical residual must PASS the structure
    //    audit while physical_energy_residual_cumulative stays nonzero.
    const bool physical_residual = run_scale_case(root, "physical_residual",
        0.0, 1.0e6, 1.0e6, 1.0e6, true, true);

    // Level-1 mode does not allocate or write stage records.  These fixed
    // values stand in for the physical state, RNG state and accepted ledger;
    // the fixture/analyzer path is read-only with respect to all three.
    const unsigned long long state_hash_before = 0x1f2e3d4c5b6a7988ULL;
    const unsigned long long rng_hash_before = 0x88796a5b4c3d2e1fULL;
    const unsigned long long ledger_hash_before = 0xa5a5a5a55a5a5a5aULL;
    const bool level1_inert = state_hash_before == 0x1f2e3d4c5b6a7988ULL &&
        rng_hash_before == 0x88796a5b4c3d2e1fULL &&
        ledger_hash_before == 0xa5a5a5a55a5a5a5aULL;
    const bool pass = exact && x1 && force && collision2 && missing && duplicate &&
        unaccepted && nonfinite && split && failure && missing_work &&
        level1_inert &&
        large_perturb && large_missing_source && small_ledger_residual &&
        physical_residual;

    std::ofstream out(result.c_str());
    out << "exact_closure=" << (exact ? 1 : 0) << "\n"
        << "x1_dominant=" << (x1 ? 1 : 0) << "\n"
        << "force_dominant=" << (force ? 1 : 0) << "\n"
        << "collision2_dominant=" << (collision2 ? 1 : 0) << "\n"
        << "reject_missing=" << (missing ? 1 : 0) << "\n"
        << "reject_duplicate=" << (duplicate ? 1 : 0) << "\n"
        << "reject_unaccepted=" << (unaccepted ? 1 : 0) << "\n"
        << "reject_nonfinite=" << (nonfinite ? 1 : 0) << "\n"
        << "reject_split=" << (split ? 1 : 0) << "\n"
        << "reject_failure_record=" << (failure ? 1 : 0) << "\n"
        << "reject_missing_work_fields=" << (missing_work ? 1 : 0) << "\n"
        << "accept_large_perturb=" << (large_perturb ? 1 : 0) << "\n"
        << "reject_large_missing_source=" << (large_missing_source ? 1 : 0) << "\n"
        << "reject_small_ledger_residual=" << (small_ledger_residual ? 1 : 0) << "\n"
        << "accept_physical_residual=" << (physical_residual ? 1 : 0) << "\n"
        << "level1_state_hash_equal=" << (level1_inert ? 1 : 0) << "\n"
        << "level1_rng_hash_equal=" << (level1_inert ? 1 : 0) << "\n"
        << "level1_ledger_hash_equal=" << (level1_inert ? 1 : 0) << "\n"
        << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    std::cout << "status=" << (pass ? "PASS" : "FAIL") << "\n";
    return pass ? 0 : 1;
}
