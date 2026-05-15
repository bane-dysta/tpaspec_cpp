// tpaspec.cpp
// C++17 single-file utility for SOS-based degenerate two-photon absorption spectra.
//
// Main workflows:
//   1) tpaspec --sos SOS.txt --out mol_tpa
//   2) tpaspec mol.fch mol.log --multiwfn /path/to/Multiwfn --out mol_tpa
//   3) tpaspec mol.fch        # auto-finds mol.log/mol.out, Multiwfn from PATH
//   4) tpaspec mol.log        # auto-finds mol.fch/mol.fchk, Multiwfn from PATH
//
// Compile:
//   g++ -std=c++17 -O2 -Wall -Wextra -pedantic tpaspec.cpp -o tpaspec
//
// Notes:
//   * The SOS.txt parser expects the Multiwfn-style format:
//       N
//       state_index energy_eV       (N lines)
//       i j mux muy muz             (dipole matrix, usually upper triangle, au)
//   * The TPA tensor uses degenerate photons, omega = E_f/2.
//   * delta_au is computed from rotational invariants of the SOS transition tensor.
//   * sigma_GM uses a normalized lineshape, so the peak GM depends on the chosen
//     HWHM/FWHM and shape. Always report these parameters.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace tpa {

constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double EH_EV = 27.211386245988;             // CODATA-ish Hartree in eV
constexpr double HC_EV_NM = 1239.8419843320026;       // h*c in eV nm
constexpr double ALPHA = 7.2973525693e-3;             // fine-structure constant
constexpr double A0_CM = 5.29177210903e-9;            // Bohr radius in cm
constexpr double C_CM_S = 2.99792458e10;              // speed of light in cm/s

using Vec3 = std::array<double, 3>;
using Cx = std::complex<double>;
using Mat3 = std::array<std::array<Cx, 3>, 3>;

std::string lower(std::string s) {
    for (char &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

bool startsWith(std::string_view s, std::string_view p) {
    return s.size() >= p.size() && s.substr(0, p.size()) == p;
}

double parseDouble(const std::string &s, const std::string &what) {
    try {
        size_t pos = 0;
        double v = std::stod(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception &) {
        throw std::runtime_error("Cannot parse " + what + " as number: " + s);
    }
}

int parseInt(const std::string &s, const std::string &what) {
    try {
        size_t pos = 0;
        int v = std::stoi(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception &) {
        throw std::runtime_error("Cannot parse " + what + " as integer: " + s);
    }
}

std::string shellQuote(const std::string &s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\"'\"'";
        else out += c;
    }
    out += "'";
    return out;
}

std::string pathNoExt(const fs::path &p) {
    fs::path q = p;
    q.replace_extension("");
    return q.string();
}

std::string csvEscape(const std::string &s) {
    bool need = s.find_first_of(",\"\n\r") != std::string::npos;
    if (!need) return s;
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

enum class Polarization { Linear, Perpendicular, Circular };
enum class Shape { Lorentzian, Gaussian };
enum class WidthKind { HWHM, FWHM };
enum class XAxis { ExcEV, PhotonEV, NM };

struct Options {
    fs::path sos;
    fs::path fch;
    fs::path log;
    std::string multiwfn = "Multiwfn";
    std::string outPrefix;

    Polarization pol = Polarization::Linear;
    Shape shape = Shape::Lorentzian;
    WidthKind widthKind = WidthKind::HWHM;
    XAxis xAxis = XAxis::NM;

    double widthEV = 0.10;      // interpreted by widthKind
    double etaSosEV = 0.0;      // imaginary damping in SOS denominators
    double nFactor = 4.0;       // common single-beam convention with S as in full SOS formula
    int points = 2000;
    int maxFinalStates = -1;    // -1 = all
    std::optional<double> gridMinEV;
    std::optional<double> gridMaxEV;

    bool dryRun = false;
    bool keepWork = false;
    bool writeGnuplot = true;
};

std::string toString(Polarization p) {
    switch (p) {
        case Polarization::Linear: return "linear";
        case Polarization::Perpendicular: return "perpendicular";
        case Polarization::Circular: return "circular";
    }
    return "unknown";
}

std::string toString(Shape s) {
    return s == Shape::Lorentzian ? "lorentzian" : "gaussian";
}

std::string toString(WidthKind w) {
    return w == WidthKind::HWHM ? "hwhm" : "fwhm";
}

std::string toString(XAxis x) {
    switch (x) {
        case XAxis::ExcEV: return "exc_ev";
        case XAxis::PhotonEV: return "photon_ev";
        case XAxis::NM: return "nm";
    }
    return "unknown";
}

Polarization parsePol(const std::string &s0) {
    std::string s = lower(s0);
    if (s == "linear" || s == "parallel" || s == "lin" || s == "para") return Polarization::Linear;
    if (s == "perpendicular" || s == "perp" || s == "orthogonal") return Polarization::Perpendicular;
    if (s == "circular" || s == "circ") return Polarization::Circular;
    throw std::runtime_error("Unknown polarization: " + s0);
}

Shape parseShape(const std::string &s0) {
    std::string s = lower(s0);
    if (s == "lorentzian" || s == "lor" || s == "l") return Shape::Lorentzian;
    if (s == "gaussian" || s == "gauss" || s == "g") return Shape::Gaussian;
    throw std::runtime_error("Unknown lineshape: " + s0);
}

WidthKind parseWidthKind(const std::string &s0) {
    std::string s = lower(s0);
    if (s == "hwhm" || s == "half") return WidthKind::HWHM;
    if (s == "fwhm" || s == "full") return WidthKind::FWHM;
    throw std::runtime_error("Unknown width kind: " + s0);
}

XAxis parseXAxis(const std::string &s0) {
    std::string s = lower(s0);
    if (s == "exc_ev" || s == "ev" || s == "energy" || s == "excitation") return XAxis::ExcEV;
    if (s == "photon_ev" || s == "photonev" || s == "photon") return XAxis::PhotonEV;
    if (s == "nm" || s == "wavelength" || s == "lambda") return XAxis::NM;
    throw std::runtime_error("Unknown x-axis: " + s0);
}

bool isSosPath(const fs::path &p) {
    std::string ext = lower(p.extension().string());
    std::string name = lower(p.filename().string());
    return ext == ".sos" || ext == ".txt" || name == "sos.txt";
}

bool isFchPath(const fs::path &p) {
    std::string ext = lower(p.extension().string());
    return ext == ".fch" || ext == ".fchk";
}

bool isLogPath(const fs::path &p) {
    std::string ext = lower(p.extension().string());
    return ext == ".log" || ext == ".out";
}

std::optional<fs::path> findSiblingWithExts(const fs::path &p, const std::vector<std::string> &exts) {
    fs::path base = p.parent_path() / p.stem();
    for (const auto &e : exts) {
        fs::path cand = base.string() + e;
        if (fs::exists(cand)) return cand;
    }
    // Also try by removing two extensions, e.g. mol.g16.log -> mol.fch
    fs::path stem2 = p.stem();
    if (!stem2.extension().empty()) {
        fs::path base2 = p.parent_path() / stem2.stem();
        for (const auto &e : exts) {
            fs::path cand = base2.string() + e;
            if (fs::exists(cand)) return cand;
        }
    }
    return std::nullopt;
}

void printUsage(std::ostream &os) {
    os << R"USAGE(tpaspec - SOS-based degenerate TPA spectrum from Multiwfn SOS.txt

Usage:
  tpaspec --sos SOS.txt [options]
  tpaspec SOS.txt [options]
  tpaspec mol.fch mol.log [options]
  tpaspec mol.fch [options]       # auto-find same-stem .log/.out
  tpaspec mol.log [options]       # auto-find same-stem .fch/.fchk

Core options:
  --sos FILE             Use existing Multiwfn SOS.txt; skip Multiwfn
  --fch FILE             Gaussian formatted checkpoint
  --log FILE             Gaussian TD log/out
  --multiwfn EXE         Multiwfn executable; default: Multiwfn from PATH
  --out PREFIX           Output prefix; default: input stem + _tpa

TPA/Spectrum options:
  --pol linear|perpendicular|circular   default: linear
  --shape lorentzian|gaussian           default: lorentzian
  --width-eV VALUE                      default: 0.10
  --width-kind hwhm|fwhm                default: hwhm
  --eta-sos-eV VALUE                    SOS denominator damping; default: 0
  --n-factor VALUE                      GM conversion N factor; default: 4
  --x nm|exc_ev|photon_ev               gnuplot x-axis; default: nm
  --points N                            curve grid points; default: 2000
  --states N                            final states to include; default: all
  --grid-min-eV VALUE                   total excitation energy grid minimum
  --grid-max-eV VALUE                   total excitation energy grid maximum

Multiwfn-run options:
  --dry-run             Show Multiwfn input script/command and exit
  --keep-work           Keep temporary Multiwfn working directory

Outputs:
  PREFIX.sticks.csv     per-state E, photon wavelength, delta_au, peak GM
  PREFIX.curve.csv      broadened curve on total excitation-energy grid
  PREFIX.gnuplot        simple gnuplot script
  PREFIX.SOS.txt        copy of the SOS.txt used/generated
  PREFIX.multiwfn.out   Multiwfn stdout, only when Multiwfn is run

Compile:
  g++ -std=c++17 -O2 -Wall -Wextra -pedantic tpaspec.cpp -o tpaspec
)USAGE";
}

Options parseArgs(int argc, char **argv) {
    Options opt;
    std::vector<fs::path> positional;

    auto needValue = [&](int &i, const std::string &name) -> std::string {
        if (i + 1 >= argc) throw std::runtime_error("Option requires a value: " + name);
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-h" || a == "--help") {
            printUsage(std::cout);
            std::exit(0);
        } else if (a == "--sos") {
            opt.sos = needValue(i, a);
        } else if (a == "--fch" || a == "--fchk") {
            opt.fch = needValue(i, a);
        } else if (a == "--log" || a == "--outlog") {
            opt.log = needValue(i, a);
        } else if (a == "--multiwfn") {
            opt.multiwfn = needValue(i, a);
        } else if (a == "--out") {
            opt.outPrefix = needValue(i, a);
        } else if (a == "--pol") {
            opt.pol = parsePol(needValue(i, a));
        } else if (a == "--shape") {
            opt.shape = parseShape(needValue(i, a));
        } else if (a == "--width-eV" || a == "--width") {
            opt.widthEV = parseDouble(needValue(i, a), a);
        } else if (a == "--width-kind") {
            opt.widthKind = parseWidthKind(needValue(i, a));
        } else if (a == "--eta-sos-eV" || a == "--eta") {
            opt.etaSosEV = parseDouble(needValue(i, a), a);
        } else if (a == "--n-factor" || a == "--N") {
            opt.nFactor = parseDouble(needValue(i, a), a);
        } else if (a == "--x" || a == "--x-axis") {
            opt.xAxis = parseXAxis(needValue(i, a));
        } else if (a == "--points") {
            opt.points = parseInt(needValue(i, a), a);
        } else if (a == "--states") {
            opt.maxFinalStates = parseInt(needValue(i, a), a);
        } else if (a == "--grid-min-eV") {
            opt.gridMinEV = parseDouble(needValue(i, a), a);
        } else if (a == "--grid-max-eV") {
            opt.gridMaxEV = parseDouble(needValue(i, a), a);
        } else if (a == "--dry-run") {
            opt.dryRun = true;
        } else if (a == "--keep-work") {
            opt.keepWork = true;
        } else if (a == "--no-gnuplot") {
            opt.writeGnuplot = false;
        } else if (startsWith(a, "--")) {
            throw std::runtime_error("Unknown option: " + a);
        } else {
            positional.emplace_back(a);
        }
    }

    for (const auto &p : positional) {
        if (isFchPath(p)) {
            if (!opt.fch.empty()) throw std::runtime_error("Multiple fch/fchk files given");
            opt.fch = p;
        } else if (isLogPath(p)) {
            if (!opt.log.empty()) throw std::runtime_error("Multiple log/out files given");
            opt.log = p;
        } else if (isSosPath(p)) {
            if (!opt.sos.empty()) throw std::runtime_error("Multiple SOS files given");
            opt.sos = p;
        } else {
            throw std::runtime_error("Cannot infer file type for positional argument: " + p.string());
        }
    }

    if (opt.widthEV <= 0.0) throw std::runtime_error("--width-eV must be positive");
    if (opt.etaSosEV < 0.0) throw std::runtime_error("--eta-sos-eV must be non-negative");
    if (opt.points < 2) throw std::runtime_error("--points must be at least 2");
    if (opt.nFactor <= 0.0) throw std::runtime_error("--n-factor must be positive");

    if (opt.sos.empty()) {
        if (!opt.fch.empty() && opt.log.empty()) {
            auto found = findSiblingWithExts(opt.fch, {".log", ".out", ".LOG", ".OUT"});
            if (found) opt.log = *found;
        }
        if (!opt.log.empty() && opt.fch.empty()) {
            auto found = findSiblingWithExts(opt.log, {".fch", ".fchk", ".FCH", ".FChk", ".FCHK"});
            if (found) opt.fch = *found;
        }
        if (opt.fch.empty() || opt.log.empty()) {
            throw std::runtime_error("Need either --sos SOS.txt, or enough fch/log input to run Multiwfn. Same-stem counterpart was not found.");
        }
    }

    if (opt.outPrefix.empty()) {
        if (!opt.sos.empty()) opt.outPrefix = pathNoExt(opt.sos) + "_tpa";
        else if (!opt.fch.empty()) opt.outPrefix = pathNoExt(opt.fch) + "_tpa";
        else opt.outPrefix = pathNoExt(opt.log) + "_tpa";
    }

    return opt;
}

struct SosData {
    int nstates = 0;
    std::vector<double> energyEV; // size nstates+1; energyEV[0]=0
    std::vector<Vec3> muFlat;     // flattened (nstates+1)*(nstates+1)

    Vec3 &mu(int i, int j) {
        return muFlat[static_cast<size_t>(i) * (nstates + 1) + static_cast<size_t>(j)];
    }
    const Vec3 &mu(int i, int j) const {
        return muFlat[static_cast<size_t>(i) * (nstates + 1) + static_cast<size_t>(j)];
    }
};

SosData readSos(const fs::path &path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open SOS file: " + path.string());

    SosData d;
    if (!(in >> d.nstates)) throw std::runtime_error("Cannot read number of states from SOS file: " + path.string());
    if (d.nstates <= 0) throw std::runtime_error("Invalid number of states in SOS file");

    d.energyEV.assign(static_cast<size_t>(d.nstates + 1), 0.0);
    d.muFlat.assign(static_cast<size_t>(d.nstates + 1) * static_cast<size_t>(d.nstates + 1), Vec3{0.0, 0.0, 0.0});

    for (int k = 0; k < d.nstates; ++k) {
        int idx = 0;
        double ev = 0.0;
        if (!(in >> idx >> ev)) throw std::runtime_error("Cannot read state energy line " + std::to_string(k + 1));
        if (idx < 1 || idx > d.nstates) throw std::runtime_error("State index out of range in energy section: " + std::to_string(idx));
        d.energyEV[idx] = ev;
    }

    int i = 0, j = 0;
    double x = 0.0, y = 0.0, z = 0.0;
    int nDip = 0;
    while (in >> i >> j >> x >> y >> z) {
        if (i < 0 || i > d.nstates || j < 0 || j > d.nstates) {
            std::ostringstream oss;
            oss << "Dipole matrix index out of range: " << i << " " << j;
            throw std::runtime_error(oss.str());
        }
        d.mu(i, j) = Vec3{x, y, z};
        d.mu(j, i) = Vec3{x, y, z};
        ++nDip;
    }
    if (!in.eof()) throw std::runtime_error("Failed while parsing dipole section in SOS file");
    if (nDip == 0) throw std::runtime_error("No dipole matrix entries found in SOS file");

    for (int k = 1; k <= d.nstates; ++k) {
        if (d.energyEV[k] <= 0.0) {
            throw std::runtime_error("Missing or non-positive excitation energy for state " + std::to_string(k));
        }
    }

    return d;
}

struct MultiwfnResult {
    fs::path sosPath;
    fs::path workDir;
    fs::path stdoutPath;
};

MultiwfnResult runMultiwfnToMakeSos(const Options &opt) {
    if (opt.fch.empty() || opt.log.empty()) {
        throw std::runtime_error("runMultiwfnToMakeSos requires both fch and log files");
    }
    if (!fs::exists(opt.fch)) throw std::runtime_error("fch/fchk file not found: " + opt.fch.string());
    if (!fs::exists(opt.log)) throw std::runtime_error("log/out file not found: " + opt.log.string());

    fs::path outPrefixPath(opt.outPrefix);
    fs::path parent = outPrefixPath.parent_path();
    if (parent.empty()) parent = fs::current_path();
    fs::create_directories(parent);

    fs::path work = fs::absolute(parent / (outPrefixPath.filename().string() + ".multiwfn_work"));
    fs::remove_all(work);
    fs::create_directories(work);

    fs::path fchCopy = work / opt.fch.filename();
    fs::path logCopy = work / opt.log.filename();
    fs::copy_file(opt.fch, fchCopy, fs::copy_options::overwrite_existing);
    fs::copy_file(opt.log, logCopy, fs::copy_options::overwrite_existing);

    // Multiwfn menu sequence used for current Multiwfn versions:
    // 18: Electronic excitation analyses
    // 5 : Calculate transition electric/magnetic dipole moments between all states
    //     then provide Gaussian TD log file
    // 3 : Generate SOS.txt
    // 0,0: return/quit
    std::ostringstream script;
    script << "18\n";
    script << "5\n";
    script << logCopy.filename().string() << "\n";
    script << "3\n";
    script << "0\n";
    script << "0\n";

    fs::path inPath = work / "multiwfn.in";
    {
        std::ofstream f(inPath);
        if (!f) throw std::runtime_error("Cannot write Multiwfn input script");
        f << script.str();
    }

    fs::path stdoutPath = work / "multiwfn.out";
    std::ostringstream cmd;
    cmd << "cd " << shellQuote(work.string()) << " && "
        << shellQuote(opt.multiwfn) << " " << shellQuote(fchCopy.filename().string())
        << " < " << shellQuote(inPath.filename().string())
        << " > " << shellQuote(stdoutPath.filename().string()) << " 2>&1";

    if (opt.dryRun) {
        std::cout << "Working directory: " << work << "\n";
        std::cout << "Multiwfn command:\n" << cmd.str() << "\n\n";
        std::cout << "Multiwfn stdin:\n" << script.str() << "\n";
        std::exit(0);
    }

    int ret = std::system(cmd.str().c_str());
    if (ret != 0) {
        std::ostringstream oss;
        oss << "Multiwfn command failed with status " << ret << ". See " << stdoutPath;
        throw std::runtime_error(oss.str());
    }

    fs::path sosGenerated = work / "SOS.txt";
    if (!fs::exists(sosGenerated)) {
        std::ostringstream oss;
        oss << "Multiwfn finished but SOS.txt was not found in " << work << ". See " << stdoutPath;
        throw std::runtime_error(oss.str());
    }

    fs::path sosOut = opt.outPrefix + std::string(".SOS.txt");
    fs::path stdoutOut = opt.outPrefix + std::string(".multiwfn.out");
    fs::copy_file(sosGenerated, sosOut, fs::copy_options::overwrite_existing);
    fs::copy_file(stdoutPath, stdoutOut, fs::copy_options::overwrite_existing);

    if (!opt.keepWork) fs::remove_all(work);

    return MultiwfnResult{sosOut, work, stdoutOut};
}

Mat3 zeroMat3() {
    Mat3 m{};
    for (auto &row : m) for (auto &v : row) v = Cx{0.0, 0.0};
    return m;
}

Mat3 tpaTensorForState(const SosData &d, int f, double etaAu) {
    Mat3 S = zeroMat3();
    const double EfAu = d.energyEV[f] / EH_EV;
    const double omega = EfAu / 2.0;
    for (int n = 1; n <= d.nstates; ++n) {
        Vec3 mu0n = d.mu(0, n);
        Vec3 munf = d.mu(n, f);

        // Difference-dipole correction for n == f:
        // <n|mu_bar|f> = <f|mu|f> - <0|mu|0>
        if (n == f) {
            Vec3 muff = d.mu(f, f);
            Vec3 mu00 = d.mu(0, 0);
            for (int a = 0; a < 3; ++a) munf[a] = muff[a] - mu00[a];
        }

        const double EnAu = d.energyEV[n] / EH_EV;
        Cx denom = Cx(EnAu - omega, -etaAu);
        if (std::abs(denom) == 0.0) {
            // Avoid exact divide-by-zero when eta=0 and a real resonance is exact.
            denom = Cx(std::numeric_limits<double>::min(), 0.0);
        }

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                double numer = mu0n[a] * munf[b] + mu0n[b] * munf[a];
                S[a][b] += numer / denom;
            }
        }
    }
    return S;
}

struct Invariants {
    double A = 0.0;
    double B = 0.0;
    double C = 0.0;
    double delta = 0.0;
};

std::tuple<double, double, double> polCoeffs(Polarization p) {
    switch (p) {
        case Polarization::Linear: return {2.0, 2.0, 2.0};
        case Polarization::Perpendicular: return {-1.0, 4.0, -1.0};
        case Polarization::Circular: return {-2.0, 3.0, 3.0};
    }
    return {2.0, 2.0, 2.0};
}

Invariants rotationalAverage(const Mat3 &S, Polarization pol) {
    Cx tr = S[0][0] + S[1][1] + S[2][2];
    Cx Bc = 0.0;
    Cx Cc = 0.0;
    for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
            Bc += S[a][b] * std::conj(S[a][b]);
            Cc += S[a][b] * std::conj(S[b][a]);
        }
    }

    auto [F, G, H] = polCoeffs(pol);
    Invariants inv;
    inv.A = std::norm(tr);
    inv.B = std::real(Bc);
    inv.C = std::real(Cc);
    inv.delta = (F * inv.A + G * inv.B + H * inv.C) / 30.0;

    // Guard against tiny negative values due only to numerical noise.
    if (inv.delta < 0.0 && inv.delta > -1e-12) inv.delta = 0.0;
    return inv;
}

struct Peak {
    int state = 0;
    double excEV = 0.0;       // total excitation energy, 2*photon energy
    double photonEV = 0.0;
    double tpaNM = 0.0;       // wavelength of one photon
    Invariants inv;
    double peakGM = 0.0;
    int nearestIntermediate = 0;
    double minDenomAbsEV = 0.0; // min |E_n - E_f/2|
};

double hwhmAu(const Options &opt) {
    const double hwhmEV = (opt.widthKind == WidthKind::HWHM) ? opt.widthEV : opt.widthEV / 2.0;
    return hwhmEV / EH_EV;
}

double etaAu(const Options &opt) {
    return opt.etaSosEV / EH_EV;
}

double lineShapeAu(double xAu, double gammaAu, Shape shape) {
    if (gammaAu <= 0.0) throw std::runtime_error("lineshape gamma must be positive");
    if (shape == Shape::Lorentzian) {
        return (1.0 / PI) * gammaAu / (xAu * xAu + gammaAu * gammaAu);
    }
    // Gaussian normalized over total excitation energy x, with HWHM=gamma.
    const double ln2 = std::log(2.0);
    return std::sqrt(ln2) / (std::sqrt(PI) * gammaAu) * std::exp(-ln2 * (xAu * xAu) / (gammaAu * gammaAu));
}

double gmFromDelta(double deltaAu, double photonEnergyAu, double lineshapeAu, double nFactor) {
    // σ[GM] = N * pi^3 * alpha * a0^5/c * omega^2 * g * delta * 1e50
    // a0 in cm, c in cm/s, g in Hartree^-1, omega in Hartree.
    const double pref = std::pow(PI, 3) * ALPHA * std::pow(A0_CM, 5) / C_CM_S * 1.0e50;
    return nFactor * pref * photonEnergyAu * photonEnergyAu * lineshapeAu * deltaAu;
}

std::vector<Peak> computePeaks(const SosData &d, const Options &opt) {
    int nfinal = d.nstates;
    if (opt.maxFinalStates > 0) nfinal = std::min(nfinal, opt.maxFinalStates);

    std::vector<Peak> peaks;
    peaks.reserve(static_cast<size_t>(nfinal));
    const double eta = etaAu(opt);
    const double gamma = hwhmAu(opt);
    const double g0 = lineShapeAu(0.0, gamma, opt.shape);

    for (int f = 1; f <= nfinal; ++f) {
        Mat3 S = tpaTensorForState(d, f, eta);
        Invariants inv = rotationalAverage(S, opt.pol);
        const double EfAu = d.energyEV[f] / EH_EV;
        const double omegaAu = EfAu / 2.0;

        Peak p;
        p.state = f;
        p.excEV = d.energyEV[f];
        p.photonEV = d.energyEV[f] / 2.0;
        p.tpaNM = HC_EV_NM / p.photonEV;
        p.inv = inv;
        p.peakGM = gmFromDelta(inv.delta, omegaAu, g0, opt.nFactor);

        p.minDenomAbsEV = std::numeric_limits<double>::infinity();
        p.nearestIntermediate = 0;
        for (int n = 1; n <= d.nstates; ++n) {
            double denom = std::abs(d.energyEV[n] - p.photonEV);
            if (denom < p.minDenomAbsEV) {
                p.minDenomAbsEV = denom;
                p.nearestIntermediate = n;
            }
        }
        peaks.push_back(p);
    }
    return peaks;
}

double spectrumAtExcEV(double excEV, const std::vector<Peak> &peaks, const Options &opt) {
    const double gamma = hwhmAu(opt);
    double total = 0.0;
    const double xAuE = excEV / EH_EV;
    for (const auto &p : peaks) {
        double EfAu = p.excEV / EH_EV;
        double omegaAu = EfAu / 2.0;
        double g = lineShapeAu(xAuE - EfAu, gamma, opt.shape);
        total += gmFromDelta(p.inv.delta, omegaAu, g, opt.nFactor);
    }
    return total;
}

std::pair<double, double> defaultGridEV(const std::vector<Peak> &peaks, const Options &opt) {
    if (peaks.empty()) throw std::runtime_error("No peaks to define grid");
    double minE = peaks.front().excEV;
    double maxE = peaks.front().excEV;
    for (const auto &p : peaks) {
        minE = std::min(minE, p.excEV);
        maxE = std::max(maxE, p.excEV);
    }
    double hwhmEV = (opt.widthKind == WidthKind::HWHM) ? opt.widthEV : opt.widthEV / 2.0;
    double pad = std::max(0.25, 8.0 * hwhmEV);
    minE = std::max(1e-6, minE - pad);
    maxE = maxE + pad;
    if (opt.gridMinEV) minE = *opt.gridMinEV;
    if (opt.gridMaxEV) maxE = *opt.gridMaxEV;
    if (!(maxE > minE)) throw std::runtime_error("Grid max must be larger than grid min");
    return {minE, maxE};
}

double xValueForPlot(double excEV, XAxis axis) {
    switch (axis) {
        case XAxis::ExcEV: return excEV;
        case XAxis::PhotonEV: return excEV / 2.0;
        case XAxis::NM: return 2.0 * HC_EV_NM / excEV;
    }
    return excEV;
}

std::string xLabel(XAxis axis) {
    switch (axis) {
        case XAxis::ExcEV: return "Excitation energy 2ω / eV";
        case XAxis::PhotonEV: return "Photon energy ω / eV";
        case XAxis::NM: return "Single-photon wavelength / nm";
    }
    return "x";
}

void writeSticksCsv(const fs::path &path, const std::vector<Peak> &peaks, const Options &opt) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write sticks CSV: " + path.string());
    out << std::setprecision(12);
    out << "state,E_exc_eV,photon_eV,tpa_nm,delta_au,A_au,B_au,C_au,peak_sigma_GM,";
    out << "nearest_intermediate_state,min_abs_En_minus_omega_eV,polarization,shape,width_eV,width_kind,eta_sos_eV,n_factor\n";
    for (const auto &p : peaks) {
        out << p.state << ',' << p.excEV << ',' << p.photonEV << ',' << p.tpaNM << ','
            << p.inv.delta << ',' << p.inv.A << ',' << p.inv.B << ',' << p.inv.C << ',' << p.peakGM << ','
            << p.nearestIntermediate << ',' << p.minDenomAbsEV << ','
            << csvEscape(toString(opt.pol)) << ',' << csvEscape(toString(opt.shape)) << ','
            << opt.widthEV << ',' << csvEscape(toString(opt.widthKind)) << ',' << opt.etaSosEV << ',' << opt.nFactor << "\n";
    }
}

void writeCurveCsv(const fs::path &path, const std::vector<Peak> &peaks, const Options &opt) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write curve CSV: " + path.string());
    out << std::setprecision(12);
    auto [minE, maxE] = defaultGridEV(peaks, opt);
    out << "x," << toString(opt.xAxis) << ",E_exc_eV,photon_eV,tpa_nm,sigma_GM\n";
    for (int i = 0; i < opt.points; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(opt.points - 1);
        double e = minE + (maxE - minE) * t;
        double x = xValueForPlot(e, opt.xAxis);
        double photon = e / 2.0;
        double nm = HC_EV_NM / photon;
        double y = spectrumAtExcEV(e, peaks, opt);
        out << x << ',' << toString(opt.xAxis) << ',' << e << ',' << photon << ',' << nm << ',' << y << "\n";
    }
}

void writeGnuplotScript(const fs::path &path, const Options &opt) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write gnuplot script: " + path.string());
    fs::path curve = fs::path(opt.outPrefix).filename().string() + ".curve.csv";
    fs::path png = fs::path(opt.outPrefix).filename().string() + ".png";
    out << "set datafile separator comma\n";
    out << "set terminal pngcairo size 1200,800\n";
    out << "set output '" << png.string() << "'\n";
    out << "set xlabel '" << xLabel(opt.xAxis) << "'\n";
    out << "set ylabel 'TPA cross section / GM'\n";
    out << "set title 'TPA spectrum: " << toString(opt.pol) << ", " << toString(opt.shape)
        << ", width=" << opt.widthEV << " eV " << toString(opt.widthKind)
        << ", N=" << opt.nFactor << "'\n";
    if (opt.xAxis == XAxis::NM) out << "set xrange [*:*] reverse\n";
    out << "plot '" << curve.string() << "' using 1:6 with lines title 'TPA'\n";
}

void writeReportJson(const fs::path &path, const fs::path &sosUsed, const SosData &d, const Options &opt) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write report JSON: " + path.string());
    out << std::setprecision(12);
    out << "{\n";
    out << "  \"sos\": \"" << sosUsed.string() << "\",\n";
    out << "  \"nstates\": " << d.nstates << ",\n";
    out << "  \"out_prefix\": \"" << opt.outPrefix << "\",\n";
    out << "  \"polarization\": \"" << toString(opt.pol) << "\",\n";
    out << "  \"shape\": \"" << toString(opt.shape) << "\",\n";
    out << "  \"width_eV\": " << opt.widthEV << ",\n";
    out << "  \"width_kind\": \"" << toString(opt.widthKind) << "\",\n";
    out << "  \"eta_sos_eV\": " << opt.etaSosEV << ",\n";
    out << "  \"n_factor\": " << opt.nFactor << ",\n";
    out << "  \"x_axis\": \"" << toString(opt.xAxis) << "\",\n";
    out << "  \"points\": " << opt.points << ",\n";
    out << "  \"hartree_eV\": " << EH_EV << ",\n";
    out << "  \"hc_eV_nm\": " << HC_EV_NM << "\n";
    out << "}\n";
}

void copySosIfNeeded(const fs::path &src, const fs::path &dst) {
    if (fs::equivalent(fs::absolute(src), fs::absolute(dst))) return;
    fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
}

} // namespace tpa

int main(int argc, char **argv) {
    try {
        using namespace tpa;
        Options opt = parseArgs(argc, argv);

        fs::path sosUsed;
        if (!opt.sos.empty()) {
            if (!fs::exists(opt.sos)) throw std::runtime_error("SOS file not found: " + opt.sos.string());
            sosUsed = opt.sos;
        } else {
            MultiwfnResult r = runMultiwfnToMakeSos(opt);
            sosUsed = r.sosPath;
        }

        SosData data = readSos(sosUsed);
        std::vector<Peak> peaks = computePeaks(data, opt);

        fs::path outPrefixPath(opt.outPrefix);
        fs::path parent = outPrefixPath.parent_path();
        if (!parent.empty()) fs::create_directories(parent);

        fs::path sosCopy = opt.outPrefix + std::string(".SOS.txt");
        if (!sosUsed.empty() && fs::exists(sosUsed)) {
            try { copySosIfNeeded(sosUsed, sosCopy); }
            catch (const std::exception &) { fs::copy_file(sosUsed, sosCopy, fs::copy_options::overwrite_existing); }
        }

        fs::path sticksPath = opt.outPrefix + std::string(".sticks.csv");
        fs::path curvePath = opt.outPrefix + std::string(".curve.csv");
        fs::path reportPath = opt.outPrefix + std::string(".report.json");
        writeSticksCsv(sticksPath, peaks, opt);
        writeCurveCsv(curvePath, peaks, opt);
        writeReportJson(reportPath, sosCopy, data, opt);

        if (opt.writeGnuplot) {
            fs::path gpPath = opt.outPrefix + std::string(".gnuplot");
            writeGnuplotScript(gpPath, opt);
        }

        double maxGM = -std::numeric_limits<double>::infinity();
        int maxState = 0;
        for (const auto &p : peaks) {
            if (p.peakGM > maxGM) { maxGM = p.peakGM; maxState = p.state; }
        }

        std::cout << "Read SOS: " << sosUsed << "\n";
        std::cout << "States: " << data.nstates << "; final states used: " << peaks.size() << "\n";
        std::cout << "Wrote: " << sticksPath << "\n";
        std::cout << "Wrote: " << curvePath << "\n";
        std::cout << "Wrote: " << reportPath << "\n";
        if (opt.writeGnuplot) std::cout << "Wrote: " << opt.outPrefix << ".gnuplot\n";
        std::cout << std::setprecision(8)
                  << "Max peak GM: state " << maxState << " = " << maxGM
                  << " GM (" << toString(opt.shape) << ", width=" << opt.widthEV
                  << " eV " << toString(opt.widthKind) << ", N=" << opt.nFactor << ")\n";

        bool warned = false;
        for (const auto &p : peaks) {
            if (p.minDenomAbsEV < std::max(opt.etaSosEV, 1e-12) || p.minDenomAbsEV < 0.05) {
                if (!warned) {
                    std::cerr << "WARNING: near-resonant SOS denominators detected; consider --eta-sos-eV.\n";
                    warned = true;
                }
                std::cerr << "  state " << p.state << ": nearest intermediate S" << p.nearestIntermediate
                          << ", |E_n - E_f/2| = " << p.minDenomAbsEV << " eV\n";
            }
        }

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n\n";
        tpa::printUsage(std::cerr);
        return 1;
    }
}
