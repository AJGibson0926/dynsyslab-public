#pragma once
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace dynsyslab::math {
// Exact integer or numerator/denominator strings. No floating-point inference
// of rationality is performed. Denominators normalize positive.
std::string canonicalRational(const std::string& value);
std::string rationalDifference(const std::string& a,const std::string& b);
std::string rationalRatio(const std::string& a,const std::string& b);
// nullopt denotes v_p(0)=+infinity. Prime arguments are checked (p<=1,000,000).
std::optional<int> primeValuation(const std::string& value,std::uint32_t prime);
std::string padicDistance(const std::string& a,const std::string& b,std::uint32_t prime);
struct PadicExpansion {
    std::uint32_t prime{};
    std::optional<int> valuation;
    // x=p^valuation * sum digits[j] p^j + O(p^(valuation+digits.size())).
    std::vector<std::uint32_t> digits;
};
PadicExpansion padicExpansion(const std::string& value,std::uint32_t prime,std::size_t precision);
bool compatiblePhaseResidues(const std::string& lower,const std::string& upper,
                            std::uint32_t prime,std::size_t lowerDepth);
// Finite prefix of the classical inverse limit of circle p-covers. Digits are
// least significant first. phase_j = base_phase/modulus + 2*pi*offsetTurns.
// Exact rational strings avoid rounding away deep branch choices.
struct PhaseCoverLevel {std::string modulus,residue,frequencyRatio,offsetTurns;};
std::vector<PhaseCoverLevel> phaseCoverTower(std::uint32_t prime,
                                            const std::vector<std::uint32_t>& digits);
struct PhaseCoverFlow {
    std::string basePhaseTurns, returnCount, prefixCarry;
    std::vector<std::uint32_t> digits;
    std::vector<PhaseCoverLevel> levels;
    std::vector<std::string> phaseTurns;
};
// Exact suspension translation at rational time, measured in base periods.
// basePhaseTurns is in [0,1); unknown higher digits are not inferred.
PhaseCoverFlow advancePhaseCover(std::uint32_t prime,const std::vector<std::uint32_t>& digits,
                                const std::string& basePhaseTurns,const std::string& elapsedBasePeriods);
struct PhaseCoverCharacter {
    std::string frequencyRatio;
    std::size_t minimumDepth{};
    // Value is exp(2*pi*i*phaseTurns); null if the prefix does not determine it.
    std::optional<std::string> phaseTurns;
};
// Characters indexed by Z[1/p]. Never select roots from unspecified digits.
PhaseCoverCharacter phaseCoverCharacter(std::uint32_t prime,
    const std::vector<std::uint32_t>& digits,const std::string& basePhaseTurns,
    const std::string& frequencyRatio);
struct PhaseCoverResonance {
    std::string baseToDriveRatio,characterRatio,driveHarmonic,branchClasses;
    std::size_t minimumDepth{};
};
// Exact kernel of r*(base/drive)+ell=0, r in Z[1/p], ell in Z.
// All solutions are integer multiples of the returned positive-r generator.
// Input is an explicitly rational positive ratio, never inferred from floats.
PhaseCoverResonance phaseCoverResonance(std::uint32_t prime,
                                      const std::string& baseToDriveRatio);
// Positive rational input. Euclidean algorithm, exact finite convergents.
struct ContinuedFraction {std::vector<std::string> terms,convergents;};
ContinuedFraction continuedFraction(const std::string& positiveRational);
struct SternBrocotRun {char direction{};std::string repetitions;};
std::vector<SternBrocotRun> sternBrocotPath(const std::string& positiveRational);
std::string rationalMediant(const std::string& a,const std::string& b);
struct FordCircle {std::string abscissa,radius;};
FordCircle fordCircle(const std::string& rational);
bool fordCirclesTangent(const std::string& a,const std::string& b);
// If a/b=A/B in reduced positive form, returns exact A*B. Its logarithm is
// the prime-weighted word metric, distinct from a p-adic metric.
std::string multiplicativeDistanceProduct(const std::string& a,const std::string& b);
}
