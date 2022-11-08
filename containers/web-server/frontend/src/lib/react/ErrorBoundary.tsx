import { ErrorOutline } from '@mui/icons-material';
import { Paper, Stack, StackProps, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { Component } from 'react';

interface ErrorBoundaryState {
  hasError: boolean;
  error: Error;
}

interface ErrorBoundaryProps {
  message: string;
  justifyErrorContent?: StackProps['justifyContent'];
}

export class ErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  constructor(props) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error) {
    // Update state so the next render will show the fallback UI.
    return { hasError: true, error };
  }
  componentDidCatch(error, errorInfo) {}
  render() {
    if (this.state.hasError) {
      return (
        <Box p={1}>
          <Stack direction="column" spacing={2}>
            <Stack direction="row" alignItems="center" justifyContent={this.props.justifyErrorContent} gap={1}>
              <ErrorOutline />
              <Typography>{this.props.message}</Typography>
            </Stack>
            {process.env.NODE_ENV === 'development' && (
              <Box maxWidth="min(100%,700px)" alignSelf="center">
                <Paper>
                  <Box p={1} maxWidth="100%" overflow="scroll">
                    <pre>
                      <code>{this.state.error?.stack.toString() ?? this.state.error.toString()}</code>
                    </pre>
                  </Box>
                </Paper>
              </Box>
            )}
          </Stack>
        </Box>
      );
    }
    return this.props.children;
  }
}
