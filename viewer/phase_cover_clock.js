// Display adapter for degree/tick controls; native phase_cover_flow owns the
// exact general rational API. Integer units avoid chart-gluing rounding errors.
function phaseCoverDisplayClock(prime, digits, initialDegrees, elapsedTicks) {
  if (![2, 3].includes(prime) || digits.length !== 2 ||
      digits.some(x => !Number.isInteger(x) || x < 0 || x >= prime) ||
      !Number.isInteger(initialDegrees) || initialDegrees < 0 || initialDegrees >= 360 ||
      !Number.isSafeInteger(elapsedTicks) || Math.abs(elapsedTicks) > 1000000) {
    throw new Error('Invalid phase-cover display clock');
  }
  const units = initialDegrees * 40 + elapsedTicks * 360;
  const turns = Math.floor(units / 14400), base = (units - turns * 14400) / 14400;
  const initial = digits[0] + prime * digits[1], size = prime * prime;
  const residue = ((initial + turns) % size + size) % size;
  const current = [residue % prime, Math.floor(residue / prime)];
  return {turns, base, residue, digits: current, prefixCarry: (initial + turns - residue) / size,
          phases: [base, (base + current[0]) / prime, (base + residue) / size]};
}
