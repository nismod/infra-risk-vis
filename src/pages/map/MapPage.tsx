import { FC } from 'react';

import { ErrorBoundary } from 'lib/react/ErrorBoundary';
import { StateEffectRoot } from 'lib/recoil/state-effects/StateEffectRoot';
import { useSyncRecoilState } from 'lib/recoil/sync-state';

import { viewState, viewStateEffect } from 'state/view';

import { MapPageDesktopLayout } from './layouts/MapPageDesktopLayout';
import { useIsMobile } from 'use-is-mobile';
import { MapPageMobileLayout } from './layouts/mobile/MapPageMobileLayout';

interface MapPageProps {
  view: string;
}

export const MapPage: FC<MapPageProps> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  const isMobile = useIsMobile();

  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      {isMobile ? <MapPageMobileLayout /> : <MapPageDesktopLayout />}
    </ErrorBoundary>
  );
};
