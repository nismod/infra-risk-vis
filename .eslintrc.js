module.exports = {
  extends: ['react-app', 'eslint:recommended'],
  rules: {
    'react/prop-types': 'off',
    'no-unused-vars': 'off',
    '@typescript-eslint/no-unused-vars': [
      'warn', 
      { 
        varsIgnorePattern: '^_',
        destructuredArrayIgnorePattern: '^_',
      }
    ],
    eqeqeq: ['warn', 'smart'],
    'react-hooks/exhaustive-deps': [
      'warn',
      {
        additionalHooks: '(useRecoilCallback|useRecoilTransaction_UNSTABLE)',
      },
    ],
  },
};
