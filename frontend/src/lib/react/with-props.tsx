import _ from 'lodash';
import { ComponentType, FC, memo } from 'react';

/**
 * Apply a set of props to a component, but still enable user to override the props later. Props are merged recursively.
 * @param ComponentClass the component to apply the props to
 * @param applyProps a set of props to apply
 * @param displayName name to assign to wrapper component for React DevTools
 * @returns a wrapper component that merges props apply at creation time and at use time
 * NOTE: need to investigate if there aren't any react issues associated with creating this sort of wrapper.
 * It sure is useful, though...
 */
export function withProps<P, AppliedT extends Partial<P>>(
  ComponentClass: ComponentType<P>,
  applyProps: AppliedT,
  displayName: string = 'WithProps',
): FC<Omit<P, keyof AppliedT> & Partial<AppliedT>> {
  const WithProps: FC<Omit<P, keyof AppliedT> & Partial<AppliedT>> = memo((props) => {
    const mergedProps = _.merge({}, applyProps, props);

    return <ComponentClass {...(mergedProps as any)} />;
  });
  WithProps.displayName = displayName;

  return WithProps;
}
