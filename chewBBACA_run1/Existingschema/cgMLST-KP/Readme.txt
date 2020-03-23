MKS Source Integrity® Browser plugin 1.0.8.1 for Total Commander
================================================================

 * License:
-----------

This software is released as freeware.


 * Disclaimer:
--------------

This software is provided "AS IS" without any warranty, either expressed or
implied, including, but not limited to, the implied warranties of
merchantability and fitness for a particular purpose. The author will not be
liable for any special, incidental, consequential or indirect damages due to
loss of data or any other reason.


 * Installation:
----------------

 1. Unzip the archive to an empty directory
 2. Choose "Configuration => Options => Operation => FS-Plugins"
 3. Click "Add"
 4. Choose mksbrowser.wfx
 5. Click "OK". You can now access the plugin in "Network Neighborhood".
 6. Start Universal Viewer
 7. Choose "Options => Configure user tools..."
 8. Click "Add..."
 9. Choose mksrevsel.exe
 10. Click "Open"
 11. Select "Copy" checkbox (See mksrevsel_config.png)
 12. Click "OK"
 13. Open "mksbrowser.ini"
 14. Set "ListerPath" to the path of Universal Viewer
 15. Save "mksbrowser.ini"
 16. An installation of MKS Source Integrity Enterprise Client is required.
 17. You have to add your MKS root projects to mksbrowser.ini. If mksbrowser.ini
     does not yet exist it is initialized during the very first plugin
     initialization.


 * Update Remarks:
------------------

 o With version 1.0.8.0 the directory of mksbrowser.ini changed from the plugin
   directory to the same directory as the initialization file (wincmd.ini) of Total
   Commander. However, if a mksbrowser.ini is found in the plugin directory this
   mksbrowser.ini is used instead.


 * Description:
---------------

MKS Source Integrity Browser is a Total Commander file system plugin to browse
MKS Source Integrity projects and to view members (by F3 or F4 on a file).
Furthermore, you can view the archive report (by ALT+ENTER on a file), view
member information (by ENTER on a file), view the project history (by
ALT+ENTER on a directory) and view the member revision difference.

If you add an MKS root project to mksbrowser.ini you need to specify the
 o HostName: IP address or machine name of the MKS Integrity Server
 o User: MKS Source Integrity user ID to connect to MKS Integrity Server with
 o Password: password for the MKS Source Integrity user ID to connect to MKS
   Integrity Server with
 o Port: TCP/IP port of the MKS Integrity Server to connect to
 o Project: name of the MKS Source Integrity root target project
 o DevPath: name of development path (if applicable)

About password encryption: If PasswordCryptKey is not an empty string the
password will be encrypted in mksbrowser.ini. You should not change the
PasswordCryptKey once you have added MKS root projects (by F7), otherwise the
decryption will fail. Be aware of that you will always be able to read the
password in plain text in mksbrowser.log if LogMKSCommands is set to 1.

You need to configure the path to the stand-alone Lister (lister.exe) in the
mksbrowser.ini file ([MKSBrowser Settings] section). The stand-alone Lister can
be downloaded from http://www.ghisler.com/lister. The Universal Viewer is
recommended for use with MKS Source Integrity Revision Selection Utility as an
external tool. The Universal Viewer can be downloaded from http://atorg.net.ru.

Further settings of the [MKSBrowser Settings] section of the mksbrowser.ini are
 o logging of MKS Source Integrity commands to mksbrowser.log
 o logging of debug and error messages to mksbrowser.log
 o usage of custom project icons
 o the MKS Integrity Server connection timeout value.

About member revision difference:
 1. Display the archive report (by ALT+ENTER on a file) in Universal Viewer
 2. Select the desired revision string
 3. Choose "Tools => MKS Source Integrity Revision Selection Utility"
 4. Choose the second revision string in the Revision Selection dialog
 5. Click "OK"


 * Limitations:
---------------

 o No add, checkin, checkout, drop, or lock operations can be performed on
   members or projects.
 o Binary members can not be viewed.
 o Members with names containing special characters can not be viewed.
 o Has only been tested with MKS Source Integrity Enterprise Edition Vs. >= 8.4.
 o Advanced search settings do not make sense.


 * ChangeLog:
-------------

 o Version 1.0.8.1 (22.10.2011)
   - fixed memory leaks
   - removed MPRESS and UPX binary compression
 o Version 1.0.8.0 (04.10.2011)
   - changed directory of mksbrowser.ini and mksbrowser.log from plugin
     directory to same directory as the initialization file (wincmd.ini)
     of Total Commander
   - added 64 bit support (for Total Commander >= 8 beta 1)
 o Version 1.0.7.0 (14.08.2007)
   - added: member revision difference
 o Version 1.0.6.1 (10.09.2006)
   - fixed: path to Lister (ListerPath) can have spaces
   - fixed: backward slashes in MKS root projects are converted to forward
     slashes when adding MKS root projects (by F7)
   - added: project settings can be saved as default when adding MKS root
     projects (by F7)
 o Version 1.0.6.0 (07.09.2006)
   - fixed assembling of the target MKS project and member names
   - fixed: any editor, lister or viewer can be used as stand-alone Lister
 o Version 1.0.5.0 (17.04.2006)
   - added dialog for adding an MKS root project to mksbrowser.ini (by F7)
   - added option PasswordCryptKey in section [MKSBrowser Settings] of
     mksbrowser.ini for password encryption / decryption
   - fixed: resorted icons in custom project icon library mksbrowser.icl
 o Version 1.0.4.0 (15.02.2006)
   - added the following command line commands:
     - config or setup: open mksbrowser.ini
     - about: display About MKS Source Integrity message box
     - servers: display message box about established MKS Server Connections
   - fixed potential access violations and buffer overflows
   - added: basic support for EF Commander 4.5
 o Version 1.0.3.6 (06.12.2005)
   - added: basic support for SpeedCommander 11: no more calls to internal
     commands of Total Commander (e.g. cm_ListInternalOnly)
   - added: set ShowAboutMKS to 0 in mksbrowser.ini by DEL or F8 on "About MKS
     Source Integrity..."
   - added: set ShowAboutMKSServerConnections to 0 in mksbrowser.ini by DEL or
     F8 on "About MKS Server Connections..."
   - fixed: do never popup About message box or message box about established
     MKS Server Connections if search is in progress
   - fixed: detection of members, subprojects and MKS messages
 o Version 1.0.3.5 Beta (31.10.2005)
   - added custom project icon library mksbrowser.icl
   - fixed: remove superfluous CR characters in member information view
 o Version 1.0.3.4 Beta (29.10.2005)
   - fixed: remove superfluous CR characters in member view
 o Version 1.0.3.3 Beta (26.10.2005)
   - fixed: skip size calculation of an MKS project / directory
   - fixed: pop error message if stand-alone Lister can not be started
   - fixed: set default port number to 7001 (was 0)
 o Version 1.0.3.2 Beta (15.10.2005)
   - added options LogErrors, ShowMessages and UseCustomProjectIcons in section
     [MKSBrowser Settings] of mksbrowser.ini
   - added: remove project section from mksbrowser.ini by DEL or F8 on an MKS
     root project directory
   - finally fixed: assembling of the target MKS project name
 o Version 1.0.3.1 Beta (14.10.2005)
   - fixed assembling of the target MKS project name
 o Version 1.0.2.0 (16.09.2005)
   - added: display message box about established MKS Server Connections
   - added: display archive report view by ALT+ENTER on a file
   - added: display project history view by ALT+ENTER on a directory
   - added: open mksbrowser.ini by ALT+ENTER on "MKS Source Integrity Browser"
     in "Network Neighborhood"
   - added section [MKSBrowser Settings] to mksbrowser.ini file with options:
     ListerPath, LogMKSCommands, MKSCommandTimeOut, ShowDisconnectMKSDlgAtExit,
     ShowAboutMKS and ShowAboutMKSServerConnections
   - changed: display member information view (ENTER on a file) in internal
     lister instead of message box
   - changed: connection timeout dialog offers possibility of continuing to
     connect
   - fixed: increased size limit for number of root projects in mksbrowser.ini
   - fixed: root project names in mksbrowser.ini must not exceed 260 characters
   - fixed: ALT+ENTER in plugin root directory lead to access violation
   - renamed mksbrowser.ini.template to mksbrowser.template.ini
 o Version 1.0.1.0 (13.09.2005)
   - added: display message boxes for first time login to MKS Integrity Server
   - added: ask for disconnecting and closing all MKS Integrity Client
     applications when exiting Total Commander
 o Version 1.0.0.1 (10.09.2005)
   - first public version


 * References:
--------------

 o Total Commander FileSystem Plugin Writer's Guide by Christian Ghisler
   - http://ghisler.fileburst.com/fsplugins/fspluginhelp2.1se.zip


 * Acknowledgments:
-------------------

 o CVSBrowser plugin by Alexandre Mayor
   - http://sourceforge.net/projects/tc-cvsbrowser
 o vicons3 - Icon set for Total Commander by StickyNomad
   - http://agiles.de/tcstuff/vicons3/vicons3_111.zip


 * Trademark and Copyright Statements:
--------------------------------------

 o MKS and MKS Source Integrity are trademarks or registered trademarks of MKS
   Inc. All other trademarks acknowledged. © 2005-2011. All rights reserved.
   - http://www.mks.com
 o Total Commander is Copyright © 1993-2011 by Christian Ghisler, Ghisler Software GmbH.
   - http://www.ghisler.com


 * Feedback:
------------

If you have problems, questions, suggestions please contact Thomas Beutlich.
 o Email: support@tbeu.de
 o URL: http://tbeu.totalcmd.net