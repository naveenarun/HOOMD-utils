function plotXML(filename)
    rotate3d on

    if ~strcmp(version('-release'), '2018a')
        warning('This script was written for Matlab 2018a; proceeding anyway');
    end

    drawSpheres = true;

    % Read xml file, check version and attributes

    raw_xml = xml2struct(filename);
    hxml = raw_xml.hoomd_xml;

    hoomd_xml_version = hxml.Attributes.version;
    if ~strcmp(hoomd_xml_version, '1.5')
        disp('This tool is intended for Hoomd XML version 1.5');
        disp(['Proceeding for version ' hoomd_xml_version '...']);
    end

    assert(isfield(hxml, 'configuration'), 'Configuration attribute missing.');

    config = hxml.configuration;

    assert(isfield(config, 'box'), 'Box missing from configuration.');
    assert(isfield(config, 'Attributes'), 'Configuration attributes missing.')

    ts = str2double(config.Attributes.time_step);

    % Load box

    box = config.box.Attributes;
    bdims = [str2double(box.lx), str2double(box.ly), str2double(box.lz)];
    tilts = [str2double(box.xy), str2double(box.xz), str2double(box.yz)];
    assert(isequal(tilts, [0 0 0]), 'Box tilts not supported.');
    
    % Load positions (required)

    assert(isfield(config, 'position'), 'Positions missing from configuration.');
    positions = config.position;
    natoms = round(str2double(positions.Attributes.num));
    coordlist_raw = positions.Text;
    coordlist = textscan(coordlist_raw, '%f %f %f', 'CollectOutput', 1);
    coordlist = coordlist{1};
    assert(size(coordlist, 1) == natoms, ['Only ' num2str(size(coordlist, 1)) ' of ' num2str(natoms) ' atoms found in XMl file.'])

    % Load types (optional)

    if ~isfield(config, 'type')
        warning('No atom type field found; assuming all atoms are same type.');
        typeColors = repmat([0],[1 natoms]);
    else
        types = textscan(config.type.Text, '%s');
        types = types{1};
        uniqueTypes = unique(types);
        typeDict = containers.Map(uniqueTypes, (1:size(uniqueTypes,1))/size(uniqueTypes,1));
        typeColors = cellfun(@(x) typeDict(x), types);
    end

    % DISPLAY/PLOT

    colormap(cool);
    mypts = scatter3(coordlist(:,1), coordlist(:,2), coordlist(:,3), [], typeColors);
    mypts.SizeData = 20;

    ptSizeCallback = @(s, e) set(mypts,'Sizedata', s.Value);
    uicontrol('Style', 'text', 'Position', [0 0 50 20], 'String', 'point size')
    SliderSize = uicontrol('Style', 'slider', 'position', [50 0 200 20], 'min', 1, 'max', 90, 'Value', 20, 'String', 'asdf', 'Callback', ptSizeCallback);
    hold on

    % Bonds

    if ~isfield(config, 'bond')
        warning('No bonds found.');
    else
        bondlist_raw = config.bond.Text;
        bondlist = textscan(bondlist_raw, '%*s %f %f', 'CollectOutput', 1);
        bondlist = bondlist{1};
        for i = 1:size(bondlist, 1)
            bondStart = coordlist(bondlist(i, 1) + 1, :);
            bondEnd = coordlist(bondlist(i, 2) + 1, :);
            bondAsArgs = mat2cell([bondStart; bondEnd]', [1, 1, 1], [2]);
            plot3(bondAsArgs{:}, 'k')
            hold on
        end
    end

    %axis([-bdims(1)/2, bdims(1)/2, -bdims(2)/5, bdims(2)/5, -bdims(3)/5, bdims(3)/5]);
    %xlim([-bdims(1)/2, bdims(1)/2]);
    axisdims = axis;
    ax_x = axisdims(2) - axisdims(1);
    ax_y = axisdims(4) - axisdims(3);
    ax_z = axisdims(6) - axisdims(5);
    ax_norm = max([ax_x, ax_y, ax_z]);
    pbaspect([ax_x, ax_y, ax_z]/ax_norm);
    title(['Timestep: ' num2str(ts)]);
    set(gca, 'FontSize', 16)
    axis(axisdims);
    axis vis3d
end

%% Utility functions from internet

function [ s ] = xml2struct( file )
%Convert xml file into a MATLAB structure
% [ s ] = xml2struct( file )
%
% A file containing:
% <XMLname attrib1="Some value">
%   <Element>Some text</Element>
%   <DifferentElement attrib2="2">Some more text</Element>
%   <DifferentElement attrib3="2" attrib4="1">Even more text</DifferentElement>
% </XMLname>
%
% Will produce:
% s.XMLname.Attributes.attrib1 = "Some value";
% s.XMLname.Element.Text = "Some text";
% s.XMLname.DifferentElement{1}.Attributes.attrib2 = "2";
% s.XMLname.DifferentElement{1}.Text = "Some more text";
% s.XMLname.DifferentElement{2}.Attributes.attrib3 = "2";
% s.XMLname.DifferentElement{2}.Attributes.attrib4 = "1";
% s.XMLname.DifferentElement{2}.Text = "Even more text";
%
% Please note that the following characters are substituted
% '-' by '_dash_', ':' by '_colon_' and '.' by '_dot_'
%
% Written by W. Falkena, ASTI, TUDelft, 21-08-2010
% Attribute parsing speed increased by 40% by A. Wanner, 14-6-2011
% Added CDATA support by I. Smirnov, 20-3-2012
%
% Modified by X. Mo, University of Wisconsin, 12-5-2012

    if (nargin < 1)
        clc;
        help xml2struct
        return
    end
    
    if isa(file, 'org.apache.xerces.dom.DeferredDocumentImpl') || isa(file, 'org.apache.xerces.dom.DeferredElementImpl')
        % input is a java xml object
        xDoc = file;
    else
        %check for existance
        if (exist(file,'file') == 0)
            %Perhaps the xml extension was omitted from the file name. Add the
            %extension and try again.
            if (isempty(strfind(file,'.xml')))
                file = [file '.xml'];
            end
            
            if (exist(file,'file') == 0)
                error(['The file ' file ' could not be found']);
            end
        end
        %read the xml file
        xDoc = xmlread(file);
    end
    
    %parse xDoc into a MATLAB structure
    s = parseChildNodes(xDoc);
    
end

% ----- Subfunction parseChildNodes -----
function [children,ptext,textflag] = parseChildNodes(theNode)
    % Recurse over node children.
    children = struct;
    ptext = struct; textflag = 'Text';
    if hasChildNodes(theNode)
        childNodes = getChildNodes(theNode);
        numChildNodes = getLength(childNodes);

        for count = 1:numChildNodes
            theChild = item(childNodes,count-1);
            [text,name,attr,childs,textflag] = getNodeData(theChild);
            
            if (~strcmp(name,'#text') && ~strcmp(name,'#comment') && ~strcmp(name,'#cdata_dash_section'))
                %XML allows the same elements to be defined multiple times,
                %put each in a different cell
                if (isfield(children,name))
                    if (~iscell(children.(name)))
                        %put existsing element into cell format
                        children.(name) = {children.(name)};
                    end
                    index = length(children.(name))+1;
                    %add new element
                    children.(name){index} = childs;
                    if(~isempty(fieldnames(text)))
                        children.(name){index} = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name){index}.('Attributes') = attr; 
                    end
                else
                    %add previously unknown (new) element to the structure
                    children.(name) = childs;
                    if(~isempty(text) && ~isempty(fieldnames(text)))
                        children.(name) = text; 
                    end
                    if(~isempty(attr)) 
                        children.(name).('Attributes') = attr; 
                    end
                end
            else
                ptextflag = 'Text';
                if (strcmp(name, '#cdata_dash_section'))
                    ptextflag = 'CDATA';
                elseif (strcmp(name, '#comment'))
                    ptextflag = 'Comment';
                end
                
                %this is the text in an element (i.e., the parentNode) 
                if (~isempty(regexprep(text.(textflag),'[\s]*','')))
                    if (~isfield(ptext,ptextflag) || isempty(ptext.(ptextflag)))
                        ptext.(ptextflag) = text.(textflag);
                    else
                        %what to do when element data is as follows:
                        %<element>Text <!--Comment--> More text</element>
                        
                        %put the text in different cells:
                        % if (~iscell(ptext)) ptext = {ptext}; end
                        % ptext{length(ptext)+1} = text;
                        
                        %just append the text
                        ptext.(ptextflag) = [ptext.(ptextflag) text.(textflag)];
                    end
                end
            end
            
        end
    end
end

% ----- Subfunction getNodeData -----
function [text,name,attr,childs,textflag] = getNodeData(theNode)
    % Create structure of node info.
    
    %make sure name is allowed as structure name
    name = toCharArray(getNodeName(theNode))';
    name = strrep(name, '-', '_dash_');
    name = strrep(name, ':', '_colon_');
    name = strrep(name, '.', '_dot_');

    attr = parseAttributes(theNode);
    if (isempty(fieldnames(attr))) 
        attr = []; 
    end
    
    %parse child nodes
    [childs,text,textflag] = parseChildNodes(theNode);
    
    if (isempty(fieldnames(childs)) && isempty(fieldnames(text)))
        %get the data of any childless nodes
        % faster than if any(strcmp(methods(theNode), 'getData'))
        % no need to try-catch (?)
        % faster than text = char(getData(theNode));
        text.(textflag) = toCharArray(getTextContent(theNode))';
    end
    
end

% ----- Subfunction parseAttributes -----
function attributes = parseAttributes(theNode)
    % Create attributes structure.

    attributes = struct;
    if hasAttributes(theNode)
       theAttributes = getAttributes(theNode);
       numAttributes = getLength(theAttributes);

       for count = 1:numAttributes
            %attrib = item(theAttributes,count-1);
            %attr_name = regexprep(char(getName(attrib)),'[-:.]','_');
            %attributes.(attr_name) = char(getValue(attrib));

            %Suggestion of Adrian Wanner
            str = toCharArray(toString(item(theAttributes,count-1)))';
            k = strfind(str,'='); 
            attr_name = str(1:(k(1)-1));
            attr_name = strrep(attr_name, '-', '_dash_');
            attr_name = strrep(attr_name, ':', '_colon_');
            attr_name = strrep(attr_name, '.', '_dot_');
            attributes.(attr_name) = str((k(1)+2):(end-1));
       end
    end
end
