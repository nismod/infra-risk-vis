import { ErrorOutline } from '@mui/icons-material';
import { Stack, StackProps, Typography } from '@mui/material';
import { Box } from '@mui/system';
import { Component } from 'react';

interface ErrorBoundaryState {
  hasError: boolean;
}

interface ErrorBoundaryProps {
  message: string;
  justifyErrorContent?: StackProps['justifyContent'];
}

export class ErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  constructor(props) {
    super(props);
    this.state = { hasError: false };
  }

  static getDerivedStateFromError(error) {
    // Update state so the next render will show the fallback UI.
    return { hasError: true };
  }
  componentDidCatch(error, errorInfo) {}
  render() {
    if (this.state.hasError) {
      return (
        <Box p={1}>
          <Stack direction="row" alignItems="center" justifyContent={this.props.justifyErrorContent} gap={1}>
            <ErrorOutline />
            <Typography>{this.props.message}</Typography>
          </Stack>
        </Box>
      );
    }
    return this.props.children;
  }
}
