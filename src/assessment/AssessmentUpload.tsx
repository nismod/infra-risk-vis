import { useDropzone } from 'react-dropzone';
import { useCallback, useMemo, useState } from 'react';
import { useRecoilState } from 'recoil';

import { Assessment } from 'config/assessment/assessment';
import { assessmentList, findIndexByID } from 'state/assessment';
import { Alert } from '@mui/material';

const baseStyle = {
  flex: 1,
  display: 'flex',
  alignItems: 'center',
  padding: '20px',
  borderWidth: 2,
  borderRadius: 2,
  borderColor: '#eeeeee',
  borderStyle: 'dashed',
  backgroundColor: '#fafafa',
  color: '#bdbdbd',
  outline: 'none',
  transition: 'border .24s ease-in-out',
};

const focusedStyle = {
  borderColor: '#2196f3',
};

const acceptStyle = {
  borderColor: '#00e676',
};

const rejectStyle = {
  borderColor: '#ff1744',
};

export function AssessmentUpload() {
  const [assessments, setAssessments] = useRecoilState(assessmentList);
  const [error, setError] = useState<string>();

  const onDrop = useCallback(
    (acceptedFiles) => {
      acceptedFiles.forEach((file) => {
        const reader = new FileReader();

        reader.onabort = () => setError('File reading was aborted.');
        reader.onerror = () => setError('File reading has failed.');
        reader.onload = () => {
          // @ts-ignore: readAsText should put string in result
          const text: string = reader.result;
          try {
            const data = JSON.parse(text);
            const newAssessment: Assessment = data.assessment;
            if (newAssessment) {
              const index = findIndexByID(assessments, newAssessment.id);
              if (index === -1) {
                setAssessments([...assessments, newAssessment]);
                setError(undefined);
              } else {
                setError(`Assessment ${newAssessment.id} already exists`);
              }
            } else {
              setError('Error reading assessment from file.');
            }
          } catch (ex) {
            setError('Error reading file.');
          }
        };
        reader.readAsText(file);
      });
    },
    [assessments],
  );

  const { getRootProps, getInputProps, isFocused, isDragAccept, isDragReject } = useDropzone({
    accept: { 'text/json': ['.json'] },
    onDrop,
  });

  const style = useMemo(
    () => ({
      ...baseStyle,
      ...(isFocused ? focusedStyle : {}),
      ...(isDragAccept ? acceptStyle : {}),
      ...(isDragReject ? rejectStyle : {}),
    }),
    [isFocused, isDragAccept, isDragReject],
  );

  return (
    <section className="container">
      <h3>Import from file</h3>
      <div {...getRootProps({ className: 'dropzone', style })}>
        <input {...getInputProps()} />
        <p style={{ margin: 0 }}>Drag an assessment file over this area, or click to select files</p>
      </div>
      {error ? (
        <Alert
          sx={{ my: 2 }}
          severity="error"
          onClose={() => {
            setError(undefined);
          }}
        >
          {error}
        </Alert>
      ) : null}
    </section>
  );
}
