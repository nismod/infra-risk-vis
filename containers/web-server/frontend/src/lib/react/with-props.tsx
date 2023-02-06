import _ from 'lodash';
import { ComponentType, FC, useMemo } from 'react';

/**
 * Apply a set of props to a component, but still enable user to override the props later. Props are merged recursively.
 * @param ComponentClass the component to apply the props to
 * @param applyProps a set of props to apply
 * @returns a wrapper component that merges props apply at creation time and at use time
 * NOTE: need to investigate if there aren't any react issues associated with creating this sort of wrapper.
 * It sure is useful, though...
 */
export function withProps<P>(ComponentClass: ComponentType<P>, applyProps: Partial<P>) {
  const WithProps: FC<P> = (props: P) => {
    const mergedProps = useMemo(() => _.merge({}, applyProps, props), [props]);
    return <ComponentClass {...mergedProps} />;
  };
  WithProps.displayName = ComponentClass.displayName;

  return WithProps;
}
