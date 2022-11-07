import { FC } from 'react';
import { useRecoilValue } from 'recoil';

import { StateEffectRoot } from '@/lib/recoil/state-effects/StateEffectRoot';

import { DamageSourceControl } from '../networks/DamageSourceControl';
import { RiskTypeControl } from './RiskTypeControl';
import { riskSourceTypeState, riskTargetTypeState, riskTargetTypeStateEffect } from './state';

export const RiskSection: FC = () => {
  const riskTarget = useRecoilValue(riskTargetTypeState);
  return (
    <>
      <StateEffectRoot state={riskTargetTypeState} effect={riskTargetTypeStateEffect} />
      <RiskTypeControl />
      {riskTarget === 'transport' && <DamageSourceControl />}
    </>
  );
};
