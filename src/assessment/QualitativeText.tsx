export const QualitativeText = ({ value }: { value: Number }) => {
  let text = 'neutral';
  const weak_threshold = 0.05;
  const strong_threshold = 0.3;
  if (value <= -strong_threshold) {
    text = 'strongly negative';
  }
  if (value >= strong_threshold) {
    text = 'strongly positive';
  }
  if (value > weak_threshold && value < strong_threshold) {
    text = 'slightly positive';
  }
  if (value < -weak_threshold && value > -strong_threshold) {
    text = 'slightly negative';
  }
  return <strong>{text}</strong>;
};
