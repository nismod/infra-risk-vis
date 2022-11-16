import { AssessmentView } from 'assessment/AssessmentView';
import { ErrorBoundary } from 'lib/react/ErrorBoundary';

export const AssessmentPage = () => (
  <ErrorBoundary message="There was a problem displaying this page.">
    <AssessmentView />
  </ErrorBoundary>
);