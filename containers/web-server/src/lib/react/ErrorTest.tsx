/**
 * This component is intended only for testing of the ErrorBoundary component.
 * @returns throws an error
 */
export function ErrorTest() {
  try {
    return null;
  } finally {
    //eslint-disable-next-line no-unsafe-finally
    throw new Error();
  }
}
