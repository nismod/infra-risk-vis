import React from 'react';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import Typography from '@material-ui/core/Typography';

import { titleCase } from './helpers';

const FeatureSidebar = (props) => {
  if (!props.feature) {
    return null
  }
  const f = props.feature.properties;

  return (
    <div className="custom-map-control top-right selected-feature">
      <Typography variant="h6">Selected Asset</Typography>
      <pre style={{display:"none"}}><code>
        {JSON.stringify(f, "", 2)}
      </code></pre>
      <List>
        {
          Object.entries(f).map(([key, value]) => (
            <ListItem>
              <ListItemText
                primary={titleCase(key.replace(/_/g, " "))}
                primaryTypographyProps={{variant:"caption"}}
                secondary={(value ===" ")? "-": value || "-"}
                />
            </ListItem>
          ))
        }
      </List>
    </div>
  )
}


export default FeatureSidebar;
