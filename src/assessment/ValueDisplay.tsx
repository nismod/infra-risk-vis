import { Slider } from "@mui/material";
import { CompactValue } from "./CompactValue";

export const ValueDisplay = ({ value }: { value: number }) => (
  <>
    <Slider
      disabled
      min={-1}
      max={1}
      marks={[
        { value: -1, label: '-' },
        { value: 0, label: '0' },
        { value: 1, label: '+' },
      ]}
      track={false}
      value={value}
    />
    <CompactValue label="Effect" value={value} />
  </>
);