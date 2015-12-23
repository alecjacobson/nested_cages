#ifndef IO_H
#define IO_H

// This file implements functions to help parse
// the command line arguments
//   - remove_all_chars_and_count
//   - legal_int
//

// Remove all characters 'c' from string 'str'.
// Also output the number of characters removed
//
// Inputs:
//   str  input string
//   c  character to remove from str
// Output:
//   removed  number of characters removed (number of appearences of 'c' in 'str')
//   str  overwritten 'str' without 'c'  
//
int remove_all_chars_and_count(
  char* str, 
  const char c);

// Check if a given string is a number
//
// Inputs:
//   str  input string
// Output:
//   bool  

bool legal_int(
  const char * str); 

#endif 