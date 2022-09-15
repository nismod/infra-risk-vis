import { TERRESTRIAL_LANDUSE_COLORS } from 'config/solutions/colors';
import { LANDUSE_HIERARCHY } from 'config/solutions/landuse-hierarchy';
import { CheckboxTree } from 'lib/controls/checkbox-tree/CheckboxTree';
import { useRecoilState } from 'recoil';
import { LayerLabel } from 'sidebar/ui/LayerLabel';
import { landuseTreeCheckboxState, landuseTreeConfig, landuseTreeExpandedState } from 'state/solutions/landuse-tree';

export const TerrestrialLandUseTree = () => {
  const [checkboxState, setCheckboxState] = useRecoilState(landuseTreeCheckboxState);
  const [expanded, setExpanded] = useRecoilState(landuseTreeExpandedState);

  return (
    <CheckboxTree
      nodes={LANDUSE_HIERARCHY}
      config={landuseTreeConfig}
      getLabel={(node) =>
        node.children ? (
          node.label
        ) : (
          <LayerLabel label={node.label} type="polygon" color={TERRESTRIAL_LANDUSE_COLORS[node.id].css} />
        )
      }
      checkboxState={checkboxState}
      onCheckboxState={setCheckboxState}
      expanded={expanded}
      onExpanded={setExpanded}
    />
  );
};
