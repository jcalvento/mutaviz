import React from 'react'
import { Header } from './header'
import { connect } from 'react-redux'
import { LOADING } from '../../reducers/appReducer'

const ConnectedHeader = ({enabled}) => {
  return (
    <Header enabled= {enabled}/>
  )
}

const mapStateToProps = (state) => {
  return {
    enabled: state.main.pageState !== LOADING
  }
}

export default connect(mapStateToProps, null)(ConnectedHeader)