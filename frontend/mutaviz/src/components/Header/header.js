import React from 'react'
import './header.css'

export const Header = ({ enabled }) => {
  return (
    <header >
      <a className= {`header_action ${enabled ? '' : 'header_action_disabled'}`} href= '/'>Mutaviz</a>
    </header>    
  )
}