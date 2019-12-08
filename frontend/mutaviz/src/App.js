import React from 'react';
import { createStore } from 'redux';
import { Provider } from 'react-redux';
import { rootReducer } from './rootReducer'
import { Mutaviz } from './components/Mutaviz'

const store = createStore(
  rootReducer,
  window.__REDUX_DEVTOOLS_EXTENSION__ &&
  window.__REDUX_DEVTOOLS_EXTENSION__()
)

const App = () => {
  return (
    <Provider store = {store}>
      <Mutaviz />
    </Provider>
  );
}

export default App;
