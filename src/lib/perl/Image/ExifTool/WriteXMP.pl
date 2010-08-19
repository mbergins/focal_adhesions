#------------------------------------------------------------------------------
# File:         WriteXMP.pl
#
# Description:  Write XMP meta information
#
# Revisions:    12/19/2004 - P. Harvey Created
#------------------------------------------------------------------------------
package Image::ExifTool::XMP;

use strict;
use Image::ExifTool qw(:DataAccess :Utils);

sub CheckXMP($$$);
sub SetPropertyPath($$;$$);
sub CaptureXMP($$$;$);

my $debug = 0;

# XMP structures (each structure is similar to a tag table so we can
# recurse through them in SetPropertyPath() as if they were tag tables)
my %xmpStruct = (
    ResourceRef => {
        NAMESPACE => 'stRef',
        InstanceID      => { },
        DocumentID      => { },
        VersionID       => { },
        RenditionClass  => { },
        RenditionParams => { },
        Manager         => { },
        ManagerVariant  => { },
        ManageTo        => { },
        ManageUI        => { },
    },
    ResourceEvent => {
        NAMESPACE => 'stEvt',
        action          => { },
        instanceID      => { },
        parameters      => { },
        softwareAgent   => { },
        when            => { },
    },
    JobRef => {
        NAMESPACE => 'stJob',
        name        => { },
        id          => { },
        url         => { },
    },
    Version => {
        NAMESPACE => 'stVer',
        comments    => { },
        event       => { Struct => 'ResourceEvent' },
        modifyDate  => { },
        modifier    => { },
        version     => { },
    },
    Thumbnail => {
        NAMESPACE => 'xapGImg',
        height      => { },
        width       => { },
       'format'     => { },
        image       => { },
    },
    IdentifierScheme => {
        NAMESPACE => 'xmpidq',
        Scheme      => { }, # qualifier for xmp:Identifier only
    },
    Dimensions => {
        NAMESPACE => 'stDim',
        w           => { },
        h           => { },
        unit        => { },
    },
    Colorant => {
        NAMESPACE => 'xapG',
        swatchName  => { },
        mode        => { },
        type        => { },
        cyan        => { },
        magenta     => { },
        yellow      => { },
        black       => { },
        red         => { },
        green       => { },
        blue        => { },
        L           => { },
        A           => { },
        B           => { },
    },
    Font => {
        NAMESPACE => 'stFnt',
        fontName    => { },
        fontFamily  => { },
        fontFace    => { },
        fontType    => { },
        versionString => { },
        composite   => { },
        fontFileName=> { },
        childFontFiles=> { List => 'Seq' },
    },
    # the following stuctures are different:  They don't have
    # their own namespaces -- instead they use the parent namespace
    Flash => {
        NAMESPACE => 'exif',
        Fired       => { },
        Return      => { },
        Mode        => { },
        Function    => { },
        RedEyeMode  => { },
    },
    OECF => {
        NAMESPACE => 'exif',
        Columns     => { },
        Rows        => { },
        Names       => { },
        Values      => { },
    },
    CFAPattern => {
        NAMESPACE => 'exif',
        Columns     => { },
        Rows        => { },
        Values      => { },
    },
    DeviceSettings => {
        NAMESPACE => 'exif',
        Columns     => { },
        Rows        => { },
        Settings    => { },
    },
    # Iptc4xmpCore structures
    ContactInfo => {
        NAMESPACE => 'Iptc4xmpCore',
        CiAdrCity   => { },
        CiAdrCtry   => { },
        CiAdrExtadr => { },
        CiAdrPcode  => { },
        CiAdrRegion => { },
        CiEmailWork => { },
        CiTelWork   => { },
        CiUrlWork   => { },
    },
    # Dynamic Media structures
    BeatSpliceStretch => {
        NAMESPACE => 'xmpDM',
        useFileBeatsMarker  => { },
        riseInDecibel       => { },
        riseInTimeDuration  => { },
    },
    Marker => {
        NAMESPACE => 'xmpDM',
        startTime   => { },
        duration    => { },
        comment     => { },
        name        => { },
        location    => { },
        target      => { },
        type        => { },
    },
    Media => {
        NAMESPACE => 'xmpDM',
        path        => { },
        track       => { },
        startTime   => { },
        duration    => { },
        managed     => { },
        webStatement=> { },
    },
    ProjectLink => {
        NAMESPACE => 'xmpDM',
        type        => { },
        path        => { },
    },
    ResampleStretch => {
        NAMESPACE => 'xmpDM',
        quality     => { },
    },
    Time => {
        NAMESPACE => 'xmpDM',
        value       => { },
        scale       => { },
    },
    Timecode => {
        NAMESPACE => 'xmpDM',
        timeValue   => { },
        timeFormat  => { },
    },
    TimeScaleStretch => {
        NAMESPACE => 'xmpDM',
        quality     => { },
        frameSize   => { },
        frameOverlappingPercentage => { },
    },
);

my $x_toolkit = "x:xmptk='Image::ExifTool $Image::ExifTool::VERSION'";
my $rdfDesc = 'rdf:Description';
#
# packet/xmp/rdf headers and trailers
#
my $pktOpen = "<?xpacket begin='\xef\xbb\xbf' id='W5M0MpCehiHzreSzNTczkc9d'?>\n";
my $xmlOpen = "<?xml version='1.0' encoding='UTF-8'?>\n";
my $xmpOpen = "<x:xmpmeta xmlns:x='$nsURI{x}' $x_toolkit>\n";
my $rdfOpen = "<rdf:RDF xmlns:rdf='$nsURI{rdf}'>\n";
my $rdfClose = "</rdf:RDF>\n";
my $xmpClose = "</x:xmpmeta>\n";
my $pktCloseW =  "<?xpacket end='w'?>"; # writable by default
my $pktCloseR =  "<?xpacket end='r'?>";

# Update XMP tag tables when this library is loaded:
# - generate all TagID's (required when writing)
# - generate PropertyPath for structure elements
# - add necessary inverse conversion routines
{
    my $mainTable = GetTagTable('Image::ExifTool::XMP::Main');
    GenerateTagIDs($mainTable);
    my ($mainTag, $ns, $tag);
    foreach $mainTag (keys %$mainTable) {
        my $mainInfo = $mainTable->{$mainTag};
        next unless ref $mainInfo eq 'HASH' and $mainInfo->{SubDirectory};
        my $table = GetTagTable($mainInfo->{SubDirectory}->{TagTable});
        $$table{WRITE_PROC} = \&WriteXMP;   # set WRITE_PROC for all tables
        GenerateTagIDs($table);
        $table->{CHECK_PROC} = \&CheckXMP;  # add our write check routine
        foreach $tag (TagTableKeys($table)) {
            my $tagInfo = $$table{$tag};
            next unless ref $tagInfo eq 'HASH';
            # must set PropertyPath now for all tags that are Struct elements
            # (normal tags will get set later if they are actually written)
            SetPropertyPath($table, $tag) if $$tagInfo{Struct};
            my $format = $$tagInfo{Writable};
        }
        # add new namespace if NAMESPACE is ns/uri pair
        next unless ref $$table{NAMESPACE};
        my $nsRef = $$table{NAMESPACE};
        # recognize as either a list or hash
        if (ref $nsRef eq 'ARRAY') {
            $ns = $$nsRef[0];
            $nsURI{$ns} = $$nsRef[1];
        } else { # must be a hash
            ($ns) = keys %$nsRef;
            $nsURI{$ns} = $$nsRef{$ns};
        }
        $$table{NAMESPACE} = $ns;
    }
}

#------------------------------------------------------------------------------
# Validate XMP packet and set read or read/write mode
# Inputs: 0) XMP data reference, 1) 'r' = read only, 'w' or undef = read/write
# Returns: true if XMP is good (and adds packet header/trailer if necessary)
sub ValidateXMP($;$)
{
    my ($xmpPt, $mode) = @_;
    unless ($$xmpPt =~ /^\0*<\0*\?\0*x\0*p\0*a\0*c\0*k\0*e\0*t/) {
        return '' unless $$xmpPt =~ /^<x(mp)?:xmpmeta/;
        # add required xpacket header/trailer
        $$xmpPt = $pktOpen . $$xmpPt . $pktCloseW;
    }
    $mode = 'w' unless $mode;
    my $end = substr($$xmpPt, -32, 32);
    # check for proper xpacket trailer and set r/w mode if necessary
    return '' unless $end =~ s/(e\0*n\0*d\0*=\0*['"]\0*)([rw])(\0*['"]\0*\?\0*>)/$1$mode$3/;
    substr($$xmpPt, -32, 32) = $end if $2 ne $mode;
    return 1;
}

#------------------------------------------------------------------------------
# Check XMP values for validity and format accordingly
# Inputs: 0) ExifTool object reference, 1) tagInfo hash reference,
#         2) raw value reference
# Returns: error string or undef (and may change value) on success
sub CheckXMP($$$)
{
    my ($exifTool, $tagInfo, $valPtr) = @_;
    # convert value from Latin if necessary
    if ($exifTool->{OPTIONS}->{Charset} eq 'Latin' and $$valPtr =~ /[\x80-\xff]/) {
        # convert from Latin to UTF-8
        my $val = Image::ExifTool::Latin2Unicode($$valPtr,'n');
        $$valPtr = Image::ExifTool::Unicode2UTF8($val,'n');
    }
    my $format = $tagInfo->{Writable};
    # if no format specified, value is a simple string
    return undef unless $format and $format ne 'string';
    if ($format eq 'rational' or $format eq 'real') {
        # make sure the value is a valid floating point number
        Image::ExifTool::IsFloat($$valPtr) or return 'Not a floating point number';
        if ($format eq 'rational') {
            $$valPtr = join('/', Image::ExifTool::Rationalize($$valPtr));
        }
    } elsif ($format eq 'integer') {
        # make sure the value is integer
        if (Image::ExifTool::IsInt($$valPtr)) {
            # no conversion required (converting to 'int' would remove leading '+')
        } elsif (Image::ExifTool::IsHex($$valPtr)) {
            $$valPtr = hex($$valPtr);
        } else {
            return 'Not an integer';
        }
    } elsif ($format eq 'date') {
        my ($y, $m, $d, $t, $tz);
        if ($$valPtr =~ /(\d{4}):(\d{2}):(\d{2}) (\d{2}:\d{2}(?::\d{2}(?:\.\d*)?)?)(.*)/) {
            ($y, $m, $d, $t, $tz) = ($1, $2, $3, $4, $5);
            $$valPtr = "$y-$m-${d}T$t";
        } elsif ($$valPtr =~ /^\s*\d{4}(:\d{2}){0,2}\s*$/) {
            # this is just a date (YYYY, YYYY-MM or YYYY-MM-DD)
            $$valPtr =~ tr/:/-/;
        } elsif ($$valPtr =~ /^\s*(\d{2}:\d{2}(?::\d{2}(?:\.\d*)?)?)(.*)\s*$/) {
            # this is just a time
            ($t, $tz) = ($1, $2);
            $$valPtr = $t;
        } else {
            return "Invalid date/time (use YYYY:MM:DD HH:MM:SS[.SS][+/-HH:MM|Z])";
        }
        if ($tz) {
            $tz =~ /^(Z|[+-]\d{2}:\d{2})$/ or return "Invalid time zone (must be Z or +/-HH:MM)";
            $$valPtr .= $tz;
        }
    } elsif ($format eq 'lang-alt') {
        # nothing to do
    } elsif ($format eq 'boolean') {
        if (not $$valPtr or $$valPtr =~ /false/i or $$valPtr =~ /^no$/i) {
            $$valPtr = 'False';
        } else {
            $$valPtr = 'True';
        }
    } elsif ($format eq '1') {
        # this is the entire XMP data block
        return 'Invalid XMP data' unless ValidateXMP($valPtr);
    } else {
        return "Unknown XMP format: $format";
    }
    return undef;   # success!
}

#------------------------------------------------------------------------------
# Get PropertyPath for specified tagInfo
# Inputs: 0) tagInfo reference
# Returns: PropertyPath string
sub GetPropertyPath($)
{
    my $tagInfo = shift;
    unless ($$tagInfo{PropertyPath}) {
        SetPropertyPath($$tagInfo{Table}, $$tagInfo{TagID});
    }
    return $$tagInfo{PropertyPath};
}

#------------------------------------------------------------------------------
# Set PropertyPath for specified tag (also for any structure elements)
# Inputs: 0) tagTable reference, 1) tagID, 2) structure reference (or undef),
#         3) property list up to this point (or undef), 4) true if tag is a list
sub SetPropertyPath($$;$$)
{
    my ($tagTablePtr, $tagID, $structPtr, $propList) = @_;
    my $table = $structPtr || $tagTablePtr;
    my $tagInfo = $$table{$tagID};
    my $ns = $$table{NAMESPACE};
    # don't override existing main table entry if already set by a Struct
    return if not $structPtr and $$tagInfo{PropertyPath};
    $ns or warn("No namespace for $tagID\n"), return;
    my (@propList, $listType);
    $propList and @propList = @$propList;
    push @propList, "$ns:$tagID";
    # lang-alt lists are handled specially, signified by Writable='lang-alt'
    if ($$tagInfo{Writable} and $$tagInfo{Writable} eq 'lang-alt') {
        $listType = 'Alt';
        # remove language code from property path if it exists
        $propList[-1] =~ s/-$$tagInfo{LangCode}$// if $$tagInfo{LangCode};
    } else {
        $listType = $$tagInfo{List};
    }
    # add required properties if this is a list
    push @propList, "rdf:$listType", 'rdf:li 000' if $listType and $listType ne '1';
    # set PropertyPath for all elements of this structure if necessary
    if ($$tagInfo{Struct}) {
        my $struct = $xmpStruct{$$tagInfo{Struct}};
        $struct or warn("No XMP $$tagInfo{Struct} structure!\n"), return;
        my $tag;
        foreach $tag (keys %$struct) {
            next if $tag eq 'NAMESPACE';
            SetPropertyPath($tagTablePtr, $tag, $struct, \@propList);
        }
    }
    # use tagInfo for combined tag name if this was a Struct
    if ($structPtr) {
        my $tagName = GetXMPTagID(\@propList);
        $$tagTablePtr{$tagName} or warn("Tag $tagName not found!\n"), return;
        $tagInfo = $$tagTablePtr{$tagName};
        # must check again for List's at this level
        if ($$tagInfo{Writable} and $$tagInfo{Writable} eq 'lang-alt') {
            $listType = 'Alt';
        } else {
            $listType = $$tagInfo{List};
        }
        push @propList, "rdf:$listType", 'rdf:li 000' if $listType and $listType ne '1';
    }
    # set property path for tagInfo in main table
    $$tagInfo{PropertyPath} = join '/', @propList;
}

#------------------------------------------------------------------------------
# Save XMP property name/value for rewriting
# Inputs: 0) ExifTool object reference
#         1) reference to array of XMP property path (last is current property)
#         2) property value, 3) optional reference to hash of property attributes
sub CaptureXMP($$$;$)
{
    my ($exifTool, $propList, $val, $attrs) = @_;
    return unless defined $val and @$propList > 2;
    if ($$propList[0] =~ /^x:x(a|m)pmeta$/ and
        $$propList[1] eq 'rdf:RDF' and
        $$propList[2] =~ /$rdfDesc( |$)/)
    {
        # no properties to save yet if this is just the description
        return unless @$propList > 3;
        # save information about this property
        my $capture = $exifTool->{XMP_CAPTURE};
        my $path = join('/', @$propList[3..$#$propList]);
        if (defined $$capture{$path}) {
            $exifTool->{XMP_ERROR} = "Duplicate XMP property: $path";
        } else {
            $$capture{$path} = [$val, $attrs || { }];
        }
    } else {
        $exifTool->{XMP_ERROR} = 'Improperly enclosed XMP property: ' . join('/',@$propList);
    }
}

#------------------------------------------------------------------------------
# Save information about resource containing blank node with nodeID
# Inputs: 0) reference to blank node information hash
#         1) reference to property list
#         2) property value
#         3) [optional] reference to attribute hash
# Notes: This routine and ProcessBlankInfo() are also used for reading information, but
#        are uncommon so are put in this file to reduce compile time for the common case
sub SaveBlankInfo($$$;$)
{
    my ($blankInfo, $propListPt, $val, $attrs) = @_;

    my $propPath = join '/', @$propListPt;
    my @ids = ($propPath =~ m{ #([^ /]*)}g);
    my $id;
    # split the property path at each nodeID
    foreach $id (@ids) {
        my ($pre, $prop, $post) = ($propPath =~ m{^(.*?)/([^/]*) #$id((/.*)?)$});
        defined $pre or warn("internal error parsing nodeID's"), next;
        # the element with the nodeID should be in the path prefix for subject
        # nodes and the path suffix for object nodes
        unless ($prop eq $rdfDesc) {
            if ($post) {
                $post = "/$prop$post";
            } else {
                $pre = "$pre/$prop";
            }
        }
        $blankInfo->{Prop}->{$id}->{Pre}->{$pre} = 1;
        if ((defined $post and length $post) or (defined $val and length $val)) {
            # save the property value and attributes for each unique path suffix
            $blankInfo->{Prop}->{$id}->{Post}->{$post} = [ $val, $attrs, $propPath ];
        }
    }
}

#------------------------------------------------------------------------------
# Process blank-node information
# Inputs: 0) ExifTool object ref, 1) tag table ref,
#         2) blank node information hash ref, 3) flag set for writing
sub ProcessBlankInfo($$$;$)
{
    my ($exifTool, $tagTablePtr, $blankInfo, $isWriting) = @_;
    $exifTool->VPrint(1, "  [Elements with nodeID set:]\n") unless $isWriting;
    my ($id, $pre, $post);
    # handle each nodeID separately
    foreach $id (sort keys %{$$blankInfo{Prop}}) {
        my $path = $blankInfo->{Prop}->{$id};
        # flag all resource names so we can warn later if some are unused
        my %unused;
        foreach $post (keys %{$path->{Post}}) {
            $unused{$post} = 1;
        }
        # combine property paths for all possible paths through this node
        foreach $pre (sort keys %{$path->{Pre}}) {
            # there will be no description for the object of a blank node
            next unless $pre =~ m{/$rdfDesc/};
            foreach $post (sort keys %{$path->{Post}}) {
                my @propList = split m{/}, "$pre$post";
                my ($val, $attrs) = @{$path->{Post}->{$post}};
                if ($isWriting) {
                    CaptureXMP($exifTool, \@propList, $val, $attrs);
                } else {
                    FoundXMP($exifTool, $tagTablePtr, \@propList, $val);
                }
                delete $unused{$post};
            }
        }
        # save information from unused properties (if RDF is malformed like f-spot output)
        if (%unused) {
            $exifTool->Options('Verbose') and $exifTool->Warn('An XMP resource is about nothing');
            foreach $post (sort keys %unused) {
                my ($val, $attrs, $propPath) = @{$path->{Post}->{$post}};
                my @propList = split m{/}, $propPath;
                if ($isWriting) {
                    CaptureXMP($exifTool, \@propList, $val, $attrs);
                } else {
                    FoundXMP($exifTool, $tagTablePtr, \@propList, $val);
                }
            }
        }
    }
}

#------------------------------------------------------------------------------
# Convert path to namespace used in file (this is a pain, but the XMP
# spec only suggests 'preferred' namespace prefixes...)
# Inputs: 0) ExifTool object reference, 1) property path
# Returns: conforming property path
sub ConformPathToNamespace($$)
{
    my ($exifTool, $path) = @_;
    my @propList = split('/',$path);
    my ($prop, $newKey);
    my $nsUsed = $exifTool->{XMP_NS};
    foreach $prop (@propList) {
        my ($ns, $tag) = $prop =~ /(.+?):(.*)/;
        next if $$nsUsed{$ns};
        my $uri = $nsURI{$ns};
        unless ($uri) {
            warn "No URI for namepace prefix $ns!\n";
            next;
        }
        my $ns2;
        foreach $ns2 (keys %$nsUsed) {
            next unless $$nsUsed{$ns2} eq $uri;
            # use the existing namespace prefix instead of ours
            $prop = "$ns2:$tag";
            last;
        }
    }
    return join('/',@propList);
}

#------------------------------------------------------------------------------
# sort tagInfo hash references by tag name
sub ByTagName
{
    return $$a{Name} cmp $$b{Name};
}

#------------------------------------------------------------------------------
# Write XMP information
# Inputs: 0) ExifTool object reference, 1) source dirInfo reference,
#         2) [optional] tag table reference
# Returns: with tag table: new XMP data (may be empty if no XMP data) or undef on error
#          without tag table: 1 on success, 0 if not valid XMP file, -1 on write error
# Notes: May set dirInfo InPlace flag to rewrite with specified DirLen
#        May set dirInfo ReadOnly flag to write as read-only XMP ('r' mode and no padding)
sub WriteXMP($$;$)
{
    my ($exifTool, $dirInfo, $tagTablePtr) = @_;
    $exifTool or return 1;    # allow dummy access to autoload this package
    my $dataPt = $$dirInfo{DataPt};
    my $dirStart = $$dirInfo{DirStart} || 0;
    my (%capture, %nsUsed, $xmpErr, $uuid);
    my $changed = 0;
    my $xmpFile = (not $tagTablePtr);   # this is an XMP data file if no $tagTablePtr
    # write XMP as preferred if this is an XMP file or a GIF file
    my $preferred = $xmpFile || ($$exifTool{FILE_TYPE} and $$exifTool{FILE_TYPE} eq 'GIF');
    my $verbose = $exifTool->Options('Verbose');
#
# extract existing XMP information into %capture hash
#
    # define hash in ExifTool object to capture XMP information (also causes
    # CaptureXMP() instead of FoundXMP() to be called from ParseXMPElement())
    #
    # The %capture hash is keyed on the complete property path beginning after
    # rdf:RDF/rdf:Description/.  The values are array references with the
    # following entries: 0) value, 1) attribute hash reference.
    $exifTool->{XMP_CAPTURE} = \%capture;
    $exifTool->{XMP_NS} = \%nsUsed;

    if ($xmpFile or ($dataPt and length $$dataPt)) {
        delete $exifTool->{XMP_ERROR};
        delete $exifTool->{XMP_UUID};
        # extract all existing XMP information (to the XMP_CAPTURE hash)
        my $success = ProcessXMP($exifTool, $dirInfo, $tagTablePtr);
        # don't continue if there is nothing to parse or if we had a parsing error
        unless ($success and not $exifTool->{XMP_ERROR}) {
            my $err = $exifTool->{XMP_ERROR} || 'Error parsing XMP';
            # may ignore this error only if we were successful
            if ($xmpFile) {
                my $raf = $$dirInfo{RAF};
                # allow empty XMP data so we can create something from nothing
                if ($success or not $raf->Seek(0,2) or $raf->Tell()) {
                    # no error message if not an XMP file
                    return 0 unless $exifTool->{XMP_ERROR};
                    if ($exifTool->Error($err, $success)) {
                        delete $exifTool->{XMP_CAPTURE};
                        return 0;
                    }
                }
            } else {
                if ($exifTool->Warn($err, $success)) {
                    delete $exifTool->{XMP_CAPTURE};
                    return undef;
                }
            }
        }
        $uuid = $exifTool->{XMP_UUID} || '';
        delete $exifTool->{XMP_ERROR};
        delete $exifTool->{XMP_UUID};
    } else {
        $uuid = '';
    }
#
# handle writing XMP as a block to XMP file
#
    if ($xmpFile) {
        my $tagInfo = $Image::ExifTool::Extra{XMP};
        if ($tagInfo and $exifTool->{NEW_VALUE}->{$tagInfo}) {
            my $rtnVal = 1;
            my $newVal = Image::ExifTool::GetNewValues($exifTool->{NEW_VALUE}->{$tagInfo});
            if (defined $newVal and length $newVal) {
                $exifTool->VPrint(0, "  Writing XMP as a block\n");
                ++$exifTool->{CHANGED};
                Write($$dirInfo{OutFile}, $newVal) or $rtnVal = -1;
            } else {
                $exifTool->Error("Can't delete all XMP from an XMP file");
            }
            delete $exifTool->{XMP_CAPTURE};
            return $rtnVal;
        }
    }
#
# delete groups in family 1 if requested
#
    if (grep /^XMP-.+$/, keys %{$exifTool->{DEL_GROUP}}) {
        my $del = $exifTool->{DEL_GROUP};
        my $path;
        foreach $path (keys %capture) {
            my @propList = split('/',$path); # get property list
            my ($tag, $namespace) = GetXMPTagID(\@propList);
            # translate namespace if necessary
            $namespace = $$xlatNamespace{$namespace} if $$xlatNamespace{$namespace};
            my $grp = "XMP-$namespace";
            my $ucg = uc $grp;
            if ($$del{$ucg} or ($$del{'XMP-*'} and not $$del{"-$ucg"})) {
                if ($verbose > 1) {
                    my $val = $capture{$path}->[0];
                    $exifTool->VPrint(1, "    - $grp:$tag = '$val'\n");
                }
                delete $capture{$path};
                ++$changed;
            }
        }
    }
#
# add, delete or change information as specified
#
    # get hash of all information we want to change
    # (sorted by tag name so alternate languages come last)
    my @tagInfoList = sort ByTagName $exifTool->GetNewTagInfoList();
    my $tagInfo;
    foreach $tagInfo (@tagInfoList) {
        next unless $exifTool->GetGroup($tagInfo, 0) eq 'XMP';
        my $tag = $tagInfo->{TagID};
        my $path = GetPropertyPath($tagInfo);
        unless ($path) {
            $exifTool->Warn("Can't write XMP:$tag (namespace unknown)");
            next;
        }
        # change our property path namespace prefixes to conform
        # to the ones used in this file
        $path = ConformPathToNamespace($exifTool, $path);
        # find existing property
        my $capList = $capture{$path};
        # MicrosoftPhoto screws up the case of some tags, so test for this
        unless ($capList) {
            my ($path2) = grep /^\Q$path\E$/i, keys %capture;
            $path2 and $capList = $capture{$path = $path2};
        }
        my $newValueHash = $exifTool->GetNewValueHash($tagInfo);
        my $overwrite = Image::ExifTool::IsOverwriting($newValueHash);
        my $writable = $$tagInfo{Writable} || '';
        my (%attrs, $deleted, $added);
        # delete existing entry if necessary
        if ($capList) {
            # take attributes from old values if they exist
            %attrs = %{$capList->[1]};
            if ($overwrite) {
                my ($delPath, @matchingPaths, $oldLang, $delLang);
                # check to see if this is an indexed list item
                if ($path =~ / /) {
                    my $pathPattern;
                    ($pathPattern = $path) =~ s/ 000/ \\d\{3\}/g;
                    @matchingPaths = sort grep(/^$pathPattern$/, keys %capture);
                } else {
                    push @matchingPaths, $path;
                }
                foreach $path (@matchingPaths) {
                    my ($val, $attrs) = @{$capture{$path}};
                    if ($overwrite < 0) {
                        # only overwrite specific values
                        next unless Image::ExifTool::IsOverwriting($newValueHash, UnescapeXML($val));
                    }
                    if ($writable eq 'lang-alt') {
                        # get original language code (lc for comparisons)
                        $oldLang = lc($$attrs{'xml:lang'} || 'x-default');
                        # delete all if deleting "x-default" or writing with no LangCode
                        # (XMP spec requires x-default language exist and be first in list)
                        if ($oldLang eq 'x-default' and not ($newValueHash->{Value} or
                            ($$tagInfo{LangCode} and $$tagInfo{LangCode} ne 'x-default')))
                        {
                            $delLang = 1;   # delete all languages
                            $overwrite = 1; # force overwrite
                        }
                        if ($$tagInfo{LangCode} and not $delLang) {
                            # only overwrite specified language
                            next unless lc($$tagInfo{LangCode}) eq $oldLang;
                        }
                    }
                    if ($verbose > 1) {
                        my $grp = $exifTool->GetGroup($tagInfo, 1);
                        my $tagName = $$tagInfo{Name};
                        $tagName =~ s/-$$tagInfo{LangCode}$// if $$tagInfo{LangCode};
                        $tagName .= '-' . $$attrs{'xml:lang'} if $$attrs{'xml:lang'};
                        $exifTool->VPrint(1, "    - $grp:$tagName = '$val'\n");
                    }
                    # save attributes and path from first deleted property
                    # so we can replace it exactly
                    unless ($delPath) {
                        %attrs = %$attrs;
                        $delPath = $path;
                    }
                    delete $capture{$path};
                    ++$changed;
                }
                next unless $delPath or $$tagInfo{List} or $oldLang;
                if ($delPath) {
                    $path = $delPath;
                    $deleted = 1;
                } else {
                    # don't change tag if we couldn't delete old copy
                    # unless this is a list or an lang-alt tag
                    next unless $$tagInfo{List} or $oldLang;
                    $path =~ m/ (\d{3})/g or warn "Internal error: no list index!\n", next;
                    $added = $1;
                }
            } elsif ($path =~ m/ (\d{3})/g) {
                $added = $1;
            }
            if (defined $added) {
                # add to end of list
                my $pos = pos($path) - 3;
                for (;;) {
                    substr($path, $pos, 3) = ++$added;
                    last unless $capture{$path};
                }
            }
        }
        # check to see if we want to create this tag
        # (create non-avoided tags in XMP data files by default)
        my $isCreating = (Image::ExifTool::IsCreating($newValueHash) or
                          ($preferred and not $$tagInfo{Avoid} and
                            not defined $$newValueHash{Shift}));

        # don't add new values unless...
            # ...tag existed before and was deleted, or we added it to a list
        next unless $deleted or defined $added or
            # ...tag didn't exist before and we are creating it
            (not $capList and $isCreating);

        # get list of new values (all done if no new values specified)
        my @newValues = Image::ExifTool::GetNewValues($newValueHash) or next;

        # set language attribute for lang-alt lists
        if ($writable eq 'lang-alt') {
            $attrs{'xml:lang'} = $$tagInfo{LangCode} || 'x-default';
            # must generate x-default entry as first entry if it didn't exist
            unless ($capList or lc($attrs{'xml:lang'}) eq 'x-default') {
                my $newValue = EscapeXML($newValues[0]);
                $capture{$path} = [ $newValue, { %attrs, 'xml:lang' => 'x-default' } ];
                if ($verbose > 1) {
                    my $tagName = $$tagInfo{Name};
                    $tagName =~ s/-$$tagInfo{LangCode}$/-x-default/;
                    my $grp = $exifTool->GetGroup($tagInfo, 1);
                    $exifTool->VPrint(1, "    + $grp:$tagName = '$newValue'\n");
                }
                $path =~ s/ 000/ 001/ or warn "Internal error: no list index!\n", next;
            }
        }
        
        # add new value(s) to %capture hash
        for (;;) {
            my $newValue = EscapeXML(shift @newValues);
            $capture{$path} = [ $newValue, \%attrs ];
            if ($verbose > 1) {
                my $grp = $exifTool->GetGroup($tagInfo, 1);
                $exifTool->VPrint(1, "    + $grp:$$tagInfo{Name} = '$newValue'\n");
            }
            ++$changed;
            last unless @newValues;
            $path =~ m/ (\d{3})/g or warn("Internal error: no list index!\n"), next;
            my $listIndex = $1;
            my $pos = pos($path) - 3;
            for (;;) {
                substr($path, $pos, 3) = ++$listIndex;
                last unless $capture{$path};
            }
            $capture{$path} and warn("Too many entries in XMP list!\n"), next;
        }
    }
    # remove the ExifTool members we created
    delete $exifTool->{XMP_CAPTURE};
    delete $exifTool->{XMP_NS};

    # return now if we didn't change anything
    unless ($changed) {
        return undef unless $xmpFile;   # just rewrite original XMP
        # get DataPt again because it may have been set by ProcessXMP
        $dataPt = $$dirInfo{DataPt};
        unless (defined $dataPt) {
            $exifTool->Error("Nothing to write");
            return 1;
        }
        Write($$dirInfo{OutFile}, $$dataPt) or return -1;
        return 1;
    }
#
# write out the new XMP information
#
    # start writing the XMP data
    my $newData = '';
    if ($$exifTool{XMP_NO_XPACKET}) {
        # write BOM if flag is set
        $newData .= "\xef\xbb\xbf" if $$exifTool{XMP_NO_XPACKET} == 2;
    } else {
        $newData .= $pktOpen;
    }
    $newData .= $xmlOpen if $$exifTool{XMP_IS_XML};
    $newData .= $xmpOpen . $rdfOpen;

    # initialize current property path list
    my @curPropList;
    my (%nsCur, $prop, $n);
    my @pathList = sort keys %capture;

    while (@pathList) {
        my $path = shift @pathList;
        my @propList = split('/',$path); # get property list
        # must open/close rdf:Description too
        unshift @propList, $rdfDesc;
        # make sure we have defined all necessary namespaces
        my (%nsNew, $newDesc);
        foreach $prop (@propList) {
            $prop =~ /(.*):/ or next;
            $1 eq 'rdf' and next;   # rdf namespace already defined
            my $nsNew = $nsUsed{$1};
            unless ($nsNew) {
                $nsNew = $nsURI{$1}; # we must have added a namespace
                unless ($nsNew) {
                    $xmpErr = "Undefined XMP namespace: $1";
                    next;
                }
            }
            $nsNew{$1} = $nsNew;
            # need a new description if any new namespaces
            $newDesc = 1 unless $nsCur{$1};
        }
        my $closeTo = 0;
        if ($newDesc) {
            # look forward to see if we will want to also open other namespaces
            # (this is necessary to keep lists from being broken if a property
            #  introduces a new namespace; plus it improves formatting)
            my ($path2, $ns2);
            foreach $path2 (@pathList) {
                my @ns2s = ($path2 =~ m{(?:^|/)([^/]+?):}g);
                my $opening = 0;
                foreach $ns2 (@ns2s) {
                    next if $ns2 eq 'rdf';
                    $nsNew{$ns2} and ++$opening, next;
                    last unless $opening and $nsURI{$ns2};
                    # also open this namespace
                    $nsNew{$ns2} = $nsURI{$ns2};
                }
                last unless $opening;
            }
        } else {
            # find first property where the current path differs from the new path
            for ($closeTo=0; $closeTo<@curPropList; ++$closeTo) {
                last unless $closeTo < @propList;
                last unless $propList[$closeTo] eq $curPropList[$closeTo];
            }
        }
        # close out properties down to the common base path
        while (@curPropList > $closeTo) {
            ($prop = pop @curPropList) =~ s/ .*//;
            $newData .= (' ' x scalar(@curPropList)) . " </$prop>\n";
        }
        if ($newDesc) {
            # open the new description
            $prop = $rdfDesc;
            %nsCur = %nsNew;            # save current namespaces
            $newData .= "\n <$prop rdf:about='$uuid'";
            foreach (sort keys %nsCur) {
                $newData .= "\n  xmlns:$_='$nsCur{$_}'";
            }
            $newData .= ">\n";
            push @curPropList, $prop;
        }
        # loop over all values for this new property
        my $capList = $capture{$path};
        my ($val, $attrs) = @$capList;
        $debug and print "$path = $val\n";
        # open new properties
        my $attr;
        for ($n=@curPropList; $n<$#propList; ++$n) {
            $prop = $propList[$n];
            push @curPropList, $prop;
            # remove list index if it exists
            $prop =~ s/ .*//;
            $attr = '';
            if ($prop ne $rdfDesc and $propList[$n+1] !~ /^rdf:/) {
                # need parseType='Resource' to avoid new 'rdf:Description'
                $attr = " rdf:parseType='Resource'";
            }
            $newData .= (' ' x scalar(@curPropList)) . "<$prop$attr>\n";
        }
        my $prop2 = pop @propList;   # get new property name
        $prop2 =~ s/ .*//;  # remove list index if it exists
        $newData .= (' ' x scalar(@curPropList)) . " <$prop2";
        # print out attributes
        foreach $attr (sort keys %$attrs) {
            my $attrVal = $$attrs{$attr};
            my $quot = ($attrVal =~ /'/) ? '"' : "'";
            $newData .= " $attr=$quot$attrVal$quot";
        }
        $newData .= ">$val</$prop2>\n";
    }
    # close off any open elements
    while ($prop = pop @curPropList) {
        $prop =~ s/ .*//;   # remove list index if it exists
        $newData .= (' ' x scalar(@curPropList)) . " </$prop>\n";
    }
#
# clean up, close out the XMP, and return our data
#
    # remove the ExifTool members we created
    delete $exifTool->{XMP_CAPTURE};
    delete $exifTool->{XMP_NS};

    $newData .= $rdfClose . $xmpClose;

    # (the XMP standard recommends writing 2k-4k of white space before the
    # packet trailer, with a newline every 100 characters)
    unless ($$exifTool{XMP_NO_XPACKET}) {
        my $pad = (' ' x 100) . "\n";
        if ($$dirInfo{InPlace}) {
            # pad to specified DirLen
            my $dirLen = $$dirInfo{DirLen} || length $$dataPt;
            my $len = length($newData) + length($pktCloseW);
            if ($len > $dirLen) {
                $exifTool->Warn('Not enough room to edit XMP in place');
                return undef;
            }
            my $num = int(($dirLen - $len) / length($pad));
            if ($num) {
                $newData .= $pad x $num;
                $len += length($pad) * $num;
            }
            $len < $dirLen and $newData .= (' ' x ($dirLen - $len - 1)) . "\n";
        } elsif (not $exifTool->Options('Compact') and
                 not $xmpFile and not $$dirInfo{ReadOnly})
        {
            $newData .= $pad x 24;
        }
        $newData .= ($$dirInfo{ReadOnly} ? $pktCloseR : $pktCloseW);
    }
    # return empty data if no properties exist
    $newData = '' unless %capture or $$dirInfo{InPlace};

    if ($xmpErr) {
        if ($xmpFile) {
            $exifTool->Error($xmpErr);
            return -1;
        }
        $exifTool->Warn($xmpErr);
        return undef;
    }
    $exifTool->{CHANGED} += $changed;
    $debug > 1 and $newData and print $newData,"\n";
    return $newData unless $xmpFile;
    Write($$dirInfo{OutFile}, $newData) or return -1;
    return 1;
}


1; # end

__END__

=head1 NAME

Image::ExifTool::WriteXMP.pl - Write XMP meta information

=head1 SYNOPSIS

These routines are autoloaded by Image::ExifTool::XMP.

=head1 DESCRIPTION

This file contains routines to write XMP metadata.

=head1 AUTHOR

Copyright 2003-2008, Phil Harvey (phil at owl.phy.queensu.ca)

This library is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=head1 SEE ALSO

L<Image::ExifTool::XMP(3pm)|Image::ExifTool::XMP>,
L<Image::ExifTool(3pm)|Image::ExifTool>

=cut
