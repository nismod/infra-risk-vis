import { Box } from '@mui/material';
import { AdaptationsSidebar } from './adaptations/AdaptationsSidebar';
import { FeatureSidebar } from './features/FeatureSidebar';
import { RegionDetails } from './regions/RegionDetails';
import { SolutionsSidebar } from './solutions/SolutionsSidebar';
import { useRecoilValue } from 'recoil';
import { showAdaptationsTableState } from 'state/adaptations';

export const DetailsSidebar = () => {
  const showAdaptationsTable = useRecoilValue(showAdaptationsTableState);
  return (
    <>
      <Box mb={2}>
        <SolutionsSidebar />
      </Box>
      <Box mb={2}>{showAdaptationsTable ? <AdaptationsSidebar /> : <FeatureSidebar />}</Box>
      <Box mb={2}>
        <RegionDetails />
      </Box>
    </>
  );
};
