import { AssessmentList } from 'assessment/AssessmentList';
import { AssessmentView } from 'assessment/AssessmentView';
import ScrollToTop from 'lib/hooks/scroll-to-top';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { useRecoilValue } from 'recoil';
import { currentAssessment } from 'state/assessment';

export const AssessmentPage = () => {
  const assessment = useRecoilValue(currentAssessment);
  return (
    <article>
      <ScrollToTop />
      <ErrorBoundary message="There was a problem displaying this page.">
        {
          assessment?
            <AssessmentView />
            : <AssessmentList />
        }
      </ErrorBoundary>
    </article>
  );
};