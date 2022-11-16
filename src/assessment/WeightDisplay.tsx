import { Slider } from "@mui/material";
import { CompactValue } from "./CompactValue";

export const WeightDisplay = ({ value, label }: { value: number, label: string }) => {
  label = label ?? "Weight";
  return (
    <>
      <Slider
        disabled
        min={0}
        max={1}
        marks={[
          { value: 0, label: '0' },
          { value: 1, label: '+' },
        ]}
        track={false}
        value={value} />
      <CompactValue label={label} value={value} />
    </>
  );
};