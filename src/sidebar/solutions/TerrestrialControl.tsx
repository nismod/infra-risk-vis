import { FormControl, FormLabel, Slider } from '@mui/material';
import {
  LandUseOption,
  LAND_USE_VALUES,
  TerrestrialLocationFilterType,
  TERRESTRIAL_LOCATION_FILTERS,
} from 'config/solutions/domains';
import { ParamChecklist } from 'lib/controls/params/ParamChecklist';
import { useRecoilState } from 'recoil';
import { InputSection } from 'sidebar/ui/InputSection';
import { terrestrialFiltersState } from 'state/solutions/terrestrial-filters';

export const TerrestrialControl = () => {
  const [terrestrialFilters, setTerrestrialFilters] = useRecoilState(terrestrialFiltersState);

  return (
    <>
      <ParamChecklist<LandUseOption>
        title="Land Use Types"
        options={LAND_USE_VALUES}
        checklistState={terrestrialFilters.landuse_desc}
        onChecklistState={(checklistState) =>
          setTerrestrialFilters({
            ...terrestrialFilters,
            landuse_desc: checklistState as Record<LandUseOption, boolean>,
          })
        }
        renderLabel={(key) => <>{key}</>}
      />
      <InputSection>
        <FormControl fullWidth>
          <FormLabel>Slope (degrees)</FormLabel>
          <Slider
            value={terrestrialFilters.slope_degrees}
            onChange={(event, value) =>
              setTerrestrialFilters({ ...terrestrialFilters, slope_degrees: value as [number, number] })
            }
            min={0}
            max={90}
            step={1}
            valueLabelDisplay="auto"
            marks
          />
        </FormControl>
      </InputSection>
      <InputSection>
        <FormControl fullWidth>
          <FormLabel>Elevation (m)</FormLabel>
          <Slider
            value={terrestrialFilters.elevation_m}
            onChange={(event, value) =>
              setTerrestrialFilters({ ...terrestrialFilters, elevation_m: value as [number, number] })
            }
            min={0}
            max={2250}
            step={10}
            valueLabelDisplay="auto"
            marks
          />
        </FormControl>
      </InputSection>
      <InputSection>
        <ParamChecklist<TerrestrialLocationFilterType>
          title="Location Attributes"
          options={[...TERRESTRIAL_LOCATION_FILTERS]}
          checklistState={terrestrialFilters.location_filters}
          onChecklistState={(checklistState) =>
            setTerrestrialFilters({ ...terrestrialFilters, location_filters: checklistState })
          }
          showAllNone={false}
          renderLabel={(key, label) => <>{label}</>}
        />
      </InputSection>
    </>
  );
};
