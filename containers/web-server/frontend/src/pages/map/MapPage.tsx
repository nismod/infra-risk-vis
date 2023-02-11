import { FC } from 'react';

import { ErrorBoundary } from '@/lib/react/ErrorBoundary';
import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';
import { useSyncRecoilState } from '@/lib/recoil/sync-state';

import { InitData } from '@/InitData';
import { viewState, viewStateEffect } from '@/state/view';
import { useIsMobile } from '@/use-is-mobile';

import { MapPageDesktopLayout } from './layouts/MapPageDesktopLayout';
import { MapPageMobileLayout } from './layouts/mobile/MapPageMobileLayout';

const MapPageLayout = () => {
  const isMobile = useIsMobile();

  return isMobile ? <MapPageMobileLayout /> : <MapPageDesktopLayout />;
};

export const MapPage: FC<{ view: string }> = ({ view }) => {
  useSyncRecoilState(viewState, view);

  return (
    <ErrorBoundary message="There was a problem displaying this page.">
      <InitData />
      <StateEffectRoot state={viewState} effect={viewStateEffect} />
      <MapPageLayout />
    </ErrorBoundary>
  );
};
