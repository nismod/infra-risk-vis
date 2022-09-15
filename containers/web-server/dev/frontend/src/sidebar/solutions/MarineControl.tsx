import { MarineLocationFilterType, MARINE_LOCATION_FILTERS } from 'config/solutions/domains';
import { ParamChecklist } from 'lib/controls/params/ParamChecklist';
import { useRecoilState } from 'recoil';
import { InputSection } from 'sidebar/ui/InputSection';
import { marineFiltersState } from 'state/solutions/marine-filters';

export const MarineControl = () => {
  const [marineFilters, setMarineFilters] = useRecoilState(marineFiltersState);

  return (
    <>
      {/* <ParamChecklist<LandUseOption>
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
      /> */}
      <InputSection>
        <ParamChecklist<MarineLocationFilterType>
          title="Find areas in proximity"
          options={[...MARINE_LOCATION_FILTERS]}
          checklistState={marineFilters.location_filters}
          onChecklistState={(checklistState) =>
            setMarineFilters({ ...marineFilters, location_filters: checklistState })
          }
          showAllNone={false}
          renderLabel={(key, label) => <>{label}</>}
        />
      </InputSection>
    </>
  );
};
